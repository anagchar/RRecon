//! Engine-independent logic: read pairs, mismatches, density, row formatting.
//!
//! Everything here is mirrored line for line in `rrecon/_pyengine.py`; the
//! column orders below must stay identical to `rrecon/schema.py`.

use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;

use crate::align::{aligned_window, CigarOp};

pub const OVERLAP_COLUMNS: &[&str] = &[
    "barcode",
    "haplotype",
    "phase_block",
    "chrom",
    "read1_name",
    "read1_start",
    "read1_end",
    "read1_strand",
    "read1_mapq",
    "read2_name",
    "read2_start",
    "read2_end",
    "read2_strand",
    "read2_mapq",
    "overlap_start",
    "overlap_end",
    "overlap_bp",
    "overlap_type",
    "seq1",
    "seq2",
    "n_mismatches",
    "mismatch_positions",
    "has_indel_in_overlap",
    "seq_match",
];
pub const SEGMENT_COLUMNS: &[&str] = &["barcode", "haplotype", "chrom", "start", "end", "depth"];
pub const BIN_COLUMNS: &[&str] = &[
    "barcode",
    "haplotype",
    "chrom",
    "bin_start",
    "bin_end",
    "n_reads",
];
pub const SUMMARY_COLUMNS: &[&str] = &[
    "barcode",
    "haplotype",
    "n_reads",
    "n_overlaps",
    "n_end_to_start",
    "n_contained",
    "n_junctions",
    "pct_junctions",
];
pub const KMER_COLUMNS: &[&str] = &["barcode", "haplotype", "junction_seq", "count"];

pub const UNTAGGED: &str = "untagged";
pub const MISSING: &str = ".";

#[derive(Clone, Debug)]
pub struct Params {
    pub barcode_tag: [u8; 2],
    pub hp_tag: Option<[u8; 2]>,
    pub ps_tag: Option<[u8; 2]>,
    pub min_mapq: u8,
    pub barcodes: Option<HashSet<String>>,
    pub bin_size: i64,
    pub junction_size: i64,
    pub overlap_tolerance: i64,
    pub max_mismatches: i64,
    /// 0 = no limit
    pub max_reads_per_group: u64,
    /// seq1/seq2/mismatch_positions are written only for overlaps up to this size
    pub max_seq_bp: i64,
    pub write_density: bool,
    pub write_overlaps: bool,
    pub write_junctions: bool,
    pub hts_threads: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            barcode_tag: *b"CB",
            hp_tag: None,
            ps_tag: None,
            min_mapq: 0,
            barcodes: None,
            bin_size: 1_000_000,
            junction_size: 9,
            overlap_tolerance: 0,
            max_mismatches: 1,
            max_reads_per_group: 0,
            max_seq_bp: 100,
            write_density: true,
            write_overlaps: true,
            write_junctions: true,
            hts_threads: 1,
        }
    }
}

impl Params {
    pub fn is_junction_size(&self, bp: i64) -> bool {
        (bp - self.junction_size).abs() <= self.overlap_tolerance
    }
    /// Sequences are always kept for junction-sized overlaps.
    pub fn seq_limit(&self) -> i64 {
        self.max_seq_bp
            .max(self.junction_size + self.overlap_tolerance)
    }
}

/// One filtered alignment.
#[derive(Clone, Debug)]
pub struct Read {
    pub name: String,
    pub start: i64,
    pub end: i64,
    pub reverse: bool,
    pub mapq: u8,
    pub ps: String,
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<u8>,
}

impl Read {
    fn key(&self) -> (i64, i64, &str) {
        (self.start, self.end, self.name.as_str())
    }
    fn strand(&self) -> char {
        if self.reverse {
            '-'
        } else {
            '+'
        }
    }
}

/// (barcode, haplotype)
pub type GroupKey = (String, String);

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GroupSummary {
    pub n_reads: u64,
    pub n_overlaps: u64,
    pub n_end_to_start: u64,
    pub n_contained: u64,
    pub n_junctions: u64,
}

impl GroupSummary {
    pub fn add(&mut self, o: &GroupSummary) {
        self.n_reads += o.n_reads;
        self.n_overlaps += o.n_overlaps;
        self.n_end_to_start += o.n_end_to_start;
        self.n_contained += o.n_contained;
        self.n_junctions += o.n_junctions;
    }
}

/// Detail of the overlap between two reads of one group.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PairDetail {
    pub seq1: Vec<u8>,
    pub seq2: Vec<u8>,
    pub n_mismatches: i64,
    pub mismatch_positions: Vec<i64>,
    pub has_indel: bool,
}

/// Overlap window, size and type of two overlapping reads, already ordered so
/// that `r1` sorts first by (start, end, name).
///
/// overlap = [r2.start, min(r1.end, r2.end)): this is also right when one read
/// lies entirely inside the other (the bug in the old script).
pub fn overlap_window(r1: &Read, r2: &Read) -> (i64, i64, bool) {
    let os = r2.start;
    let oe = r1.end.min(r2.end);
    let contained = r1.start == r2.start || r2.end <= r1.end;
    (os, oe, contained)
}

pub fn pair_detail(r1: &Read, r2: &Read, os: i64, oe: i64) -> PairDetail {
    let (s1, ins1) = aligned_window(r1.start, &r1.cigar, &r1.seq, os, oe);
    let (s2, ins2) = aligned_window(r2.start, &r2.cigar, &r2.seq, os, oe);
    let mut positions = Vec::new();
    for (k, (a, b)) in s1.iter().zip(s2.iter()).enumerate() {
        if a != b {
            positions.push(k as i64);
        }
    }
    let has_indel = ins1 || ins2 || s1.contains(&b'-') || s2.contains(&b'-');
    PairDetail {
        n_mismatches: positions.len() as i64,
        seq1: s1,
        seq2: s2,
        mismatch_positions: positions,
        has_indel,
    }
}

/// Per-chromosome accumulator for pairs: writes overlap and junction rows and
/// counts per group.
pub struct PairSink<'a> {
    pub params: &'a Params,
    pub chrom: &'a str,
    pub overlaps: Option<&'a mut dyn std::io::Write>,
    pub junctions: Option<&'a mut dyn std::io::Write>,
    pub groups: &'a mut BTreeMap<GroupKey, GroupSummary>,
    pub kmers: &'a mut BTreeMap<(String, String, String), u64>,
    line: String,
}

impl<'a> PairSink<'a> {
    pub fn new(
        params: &'a Params,
        chrom: &'a str,
        overlaps: Option<&'a mut dyn std::io::Write>,
        junctions: Option<&'a mut dyn std::io::Write>,
        groups: &'a mut BTreeMap<GroupKey, GroupSummary>,
        kmers: &'a mut BTreeMap<(String, String, String), u64>,
    ) -> Self {
        Self {
            params,
            chrom,
            overlaps,
            junctions,
            groups,
            kmers,
            line: String::with_capacity(512),
        }
    }

    /// Record the overlap between two reads of `group` (in any order).
    pub fn pair(&mut self, group: &GroupKey, a: &Read, b: &Read) -> std::io::Result<()> {
        let (r1, r2) = if a.key() <= b.key() { (a, b) } else { (b, a) };
        let (os, oe, contained) = overlap_window(r1, r2);
        let bp = oe - os;
        let p = self.params;
        let jsize = p.is_junction_size(bp);
        let need_detail = self.overlaps.is_some() || (jsize && self.junctions.is_some());
        let detail = if need_detail {
            Some(pair_detail(r1, r2, os, oe))
        } else {
            None
        };
        let is_junction = jsize
            && detail
                .as_ref()
                .is_some_and(|d| d.n_mismatches <= p.max_mismatches);

        let g = self.groups.entry(group.clone()).or_default();
        g.n_overlaps += 1;
        if contained {
            g.n_contained += 1;
        } else {
            g.n_end_to_start += 1;
        }
        if is_junction {
            g.n_junctions += 1;
        }

        let Some(d) = detail else { return Ok(()) };
        self.line.clear();
        format_pair(&mut self.line, p, group, self.chrom, r1, r2, os, oe, contained, &d);
        if let Some(w) = self.overlaps.as_mut() {
            w.write_all(self.line.as_bytes())?;
        }
        if is_junction {
            if let Some(w) = self.junctions.as_mut() {
                w.write_all(self.line.as_bytes())?;
            }
            let kmer = String::from_utf8_lossy(&d.seq1).into_owned();
            *self
                .kmers
                .entry((group.0.clone(), group.1.clone(), kmer))
                .or_insert(0) += 1;
        }
        Ok(())
    }
}

pub fn phase_block(ps1: &str, ps2: &str) -> String {
    if ps1 == ps2 {
        ps1.to_string()
    } else {
        format!("{ps1}/{ps2}")
    }
}

#[allow(clippy::too_many_arguments)]
fn format_pair(
    out: &mut String,
    p: &Params,
    group: &GroupKey,
    chrom: &str,
    r1: &Read,
    r2: &Read,
    os: i64,
    oe: i64,
    contained: bool,
    d: &PairDetail,
) {
    let bp = oe - os;
    let show_seq = bp <= p.seq_limit();
    let seq_match = d.n_mismatches <= p.max_mismatches && !d.has_indel;
    let _ = write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
        group.0,
        group.1,
        phase_block(&r1.ps, &r2.ps),
        chrom,
        r1.name,
        r1.start,
        r1.end,
        r1.strand(),
        r1.mapq,
        r2.name,
        r2.start,
        r2.end,
        r2.strand(),
        r2.mapq,
        os,
        oe,
        bp,
        if contained { "contained" } else { "end_to_start" },
    );
    if show_seq {
        out.push_str(std::str::from_utf8(&d.seq1).unwrap_or(MISSING));
        out.push('\t');
        out.push_str(std::str::from_utf8(&d.seq2).unwrap_or(MISSING));
    } else {
        out.push_str(".\t.");
    }
    let _ = write!(out, "\t{}\t", d.n_mismatches);
    if show_seq && !d.mismatch_positions.is_empty() {
        for (i, k) in d.mismatch_positions.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "{k}");
        }
    } else {
        out.push_str(MISSING);
    }
    let _ = writeln!(out, "\t{}\t{}", d.has_indel, seq_match);
}

/// Depth segments (depth >= 1) of a set of half-open intervals; adjacent
/// stretches with the same depth are merged.
pub fn depth_segments(intervals: &[(i64, i64)]) -> Vec<(i64, i64, i64)> {
    let mut events: Vec<(i64, i64)> = Vec::with_capacity(intervals.len() * 2);
    for &(s, e) in intervals {
        events.push((s, 1));
        events.push((e, -1));
    }
    events.sort_unstable();
    // net change per position, dropping positions where nothing changes
    let mut net: Vec<(i64, i64)> = Vec::new();
    for (pos, d) in events {
        match net.last_mut() {
            Some(last) if last.0 == pos => last.1 += d,
            _ => net.push((pos, d)),
        }
    }
    net.retain(|&(_, d)| d != 0);
    let mut out = Vec::new();
    let mut depth = 0;
    for w in 0..net.len() {
        depth += net[w].1;
        if depth > 0 && w + 1 < net.len() {
            out.push((net[w].0, net[w + 1].0, depth));
        }
    }
    out
}

pub fn write_segments(
    w: &mut dyn std::io::Write,
    group: &GroupKey,
    chrom: &str,
    intervals: &[(i64, i64)],
) -> std::io::Result<()> {
    for (s, e, d) in depth_segments(intervals) {
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", group.0, group.1, chrom, s, e, d)?;
    }
    Ok(())
}

pub fn write_bins(
    w: &mut dyn std::io::Write,
    group: &GroupKey,
    chrom: &str,
    chrom_len: i64,
    bin_size: i64,
    counts: &BTreeMap<i64, u64>,
) -> std::io::Result<()> {
    for (&b, &n) in counts {
        let s = b * bin_size;
        let e = ((b + 1) * bin_size).min(chrom_len.max(s + 1));
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", group.0, group.1, chrom, s, e, n)?;
    }
    Ok(())
}

pub fn format_pct(n: u64, d: u64) -> String {
    if d == 0 {
        "0.0000".into()
    } else {
        format!("{:.4}", 100.0 * n as f64 / d as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::OpKind::*;

    fn read(name: &str, start: i64, seq: &str) -> Read {
        Read {
            name: name.into(),
            start,
            end: start + seq.len() as i64,
            reverse: false,
            mapq: 60,
            ps: MISSING.into(),
            cigar: vec![CigarOp::new(Match, seq.len() as u32)],
            seq: seq.as_bytes().to_vec(),
        }
    }

    #[test]
    fn end_to_start_window() {
        let a = read("a", 100, "AAAAAAAAAAAAAAAAAAAA"); // 100..120
        let b = read("b", 111, "AAAAAAAAAAAAAAAAAAAA"); // 111..131
        let (os, oe, c) = overlap_window(&a, &b);
        assert_eq!((os, oe, c), (111, 120, false));
    }

    #[test]
    fn contained_read_uses_inner_end() {
        // the old script's bug: inner read must give min(ends) - inner start
        let a = read("a", 100, &"A".repeat(100)); // 100..200
        let b = read("b", 120, &"A".repeat(30)); // 120..150
        let (os, oe, c) = overlap_window(&a, &b);
        assert_eq!((os, oe, c), (120, 150, true));
    }

    #[test]
    fn same_start_is_contained() {
        let a = read("a", 100, &"A".repeat(10));
        let b = read("b", 100, &"A".repeat(30));
        assert_eq!(overlap_window(&a, &b), (100, 110, true));
    }

    #[test]
    fn mismatches_are_counted_by_offset() {
        let a = read("a", 0, "ACGTACGTA");
        let b = read("b", 0, "ACCTACGTT");
        let d = pair_detail(&a, &b, 0, 9);
        assert_eq!(d.n_mismatches, 2);
        assert_eq!(d.mismatch_positions, vec![2, 8]);
        assert!(!d.has_indel);
    }

    #[test]
    fn deletion_counts_as_indel_and_mismatch() {
        let a = read("a", 0, "ACGTACGTA");
        let mut b = read("b", 0, "ACGTCGTA");
        b.cigar = vec![
            CigarOp::new(Match, 4),
            CigarOp::new(Del, 1),
            CigarOp::new(Match, 4),
        ];
        b.end = 9;
        let d = pair_detail(&a, &b, 0, 9);
        assert!(d.has_indel);
        assert_eq!(d.mismatch_positions, vec![4]);
    }

    #[test]
    fn segments_merge_and_split() {
        let segs = depth_segments(&[(100, 300), (150, 350), (200, 250)]);
        assert_eq!(
            segs,
            vec![
                (100, 150, 1),
                (150, 200, 2),
                (200, 250, 3),
                (250, 300, 2),
                (300, 350, 1)
            ]
        );
        // touching reads give one depth-1 segment
        assert_eq!(depth_segments(&[(0, 10), (10, 20)]), vec![(0, 20, 1)]);
    }

    #[test]
    fn pct_format() {
        assert_eq!(format_pct(1, 3), "33.3333");
        assert_eq!(format_pct(0, 0), "0.0000");
    }

    #[test]
    fn phase_blocks() {
        assert_eq!(phase_block("100", "100"), "100");
        assert_eq!(phase_block("100", "."), "100/.");
    }
}
