//! Reading the BAM: filters, grouping, the per-chromosome sweep and the
//! output files. Chromosomes run in parallel on a rayon pool; each writes
//! plain-text part files that are joined, in header order, into one gzip
//! stream per output at the end.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read as _, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use flate2::write::GzEncoder;
use flate2::Compression;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::alignment::record::data::field::Value;
use rayon::prelude::*;

use crate::align::{reference_span, CigarOp, OpKind};
use crate::core::*;

/// Read filtering counts, in the order they are applied.
#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub total: u64,
    pub unmapped: u64,
    pub secondary_or_supplementary: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub missing_barcode: u64,
    pub other_barcode: u64,
    pub no_alignment_span: u64,
    pub kept: u64,
}

impl Stats {
    pub fn add(&mut self, o: &Stats) {
        self.total += o.total;
        self.unmapped += o.unmapped;
        self.secondary_or_supplementary += o.secondary_or_supplementary;
        self.duplicate += o.duplicate;
        self.low_mapq += o.low_mapq;
        self.missing_barcode += o.missing_barcode;
        self.other_barcode += o.other_barcode;
        self.no_alignment_span += o.no_alignment_span;
        self.kept += o.kept;
    }
    pub fn to_vec(&self) -> Vec<(String, u64)> {
        vec![
            ("total".into(), self.total),
            ("unmapped".into(), self.unmapped),
            ("secondary_or_supplementary".into(), self.secondary_or_supplementary),
            ("duplicate".into(), self.duplicate),
            ("low_mapq".into(), self.low_mapq),
            ("missing_barcode".into(), self.missing_barcode),
            ("other_barcode".into(), self.other_barcode),
            ("no_alignment_span".into(), self.no_alignment_span),
            ("kept".into(), self.kept),
        ]
    }
}

#[derive(Default)]
pub struct ChromOut {
    pub stats: Stats,
    pub groups: BTreeMap<GroupKey, GroupSummary>,
    pub kmers: BTreeMap<(String, String, String), u64>,
    /// (barcode, haplotype, n_reads) of groups over --max-reads-per-group
    pub skipped: Vec<(String, String, u64)>,
}

pub struct RunOut {
    pub stats: Stats,
    pub skipped: Vec<(String, String, String, u64)>,
    pub files: Vec<(String, PathBuf)>,
}

enum Index {
    Bai(bam::bai::Index),
    Csi(noodles_csi::Index),
}

fn load_index(bam_path: &Path) -> io::Result<Index> {
    let s = bam_path.as_os_str().to_string_lossy();
    let bai = PathBuf::from(format!("{s}.bai"));
    if bai.exists() {
        return bam::bai::fs::read(bai).map(Index::Bai);
    }
    let csi = PathBuf::from(format!("{s}.csi"));
    if csi.exists() {
        return noodles_csi::fs::read(csi).map(Index::Csi);
    }
    let alt = bam_path.with_extension("bai");
    if alt.exists() {
        return bam::bai::fs::read(alt).map(Index::Bai);
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("no .bai or .csi index next to {s}; run `samtools index` first"),
    ))
}

fn value_to_string(v: &Value) -> Option<String> {
    Some(match v {
        Value::Character(c) => (*c as char).to_string(),
        Value::Int8(n) => n.to_string(),
        Value::UInt8(n) => n.to_string(),
        Value::Int16(n) => n.to_string(),
        Value::UInt16(n) => n.to_string(),
        Value::Int32(n) => n.to_string(),
        Value::UInt32(n) => n.to_string(),
        Value::Float(f) => f.to_string(),
        Value::String(s) | Value::Hex(s) => String::from_utf8_lossy(s).into_owned(),
        Value::Array(_) => return None,
    })
}

fn get_tag(record: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    match record.data().get(tag) {
        Some(Ok(v)) => value_to_string(&v),
        _ => None,
    }
}

fn haplotype(record: &bam::Record, tag: Option<&[u8; 2]>) -> String {
    match tag.and_then(|t| get_tag(record, t)) {
        Some(v) if v == "1" || v == "2" => v,
        _ => UNTAGGED.to_string(),
    }
}

fn cigar_ops(record: &bam::Record) -> io::Result<Vec<CigarOp>> {
    record
        .cigar()
        .iter()
        .map(|r| {
            r.map(|op| {
                use noodles_sam::alignment::record::cigar::op::Kind;
                let kind = match op.kind() {
                    Kind::Match => OpKind::Match,
                    Kind::Insertion => OpKind::Ins,
                    Kind::Deletion => OpKind::Del,
                    Kind::Skip => OpKind::Skip,
                    Kind::SoftClip => OpKind::Soft,
                    Kind::HardClip => OpKind::Hard,
                    Kind::Pad => OpKind::Pad,
                    Kind::SequenceMatch => OpKind::SeqMatch,
                    Kind::SequenceMismatch => OpKind::SeqMismatch,
                };
                CigarOp::new(kind, op.len() as u32)
            })
        })
        .collect()
}

/// Apply the read filter. Returns the read's group and the read itself
/// (with sequence only when `with_seq`), or None when it is dropped.
fn extract(
    record: &bam::Record,
    p: &Params,
    stats: &mut Stats,
    with_seq: bool,
) -> io::Result<Option<(GroupKey, Read)>> {
    stats.total += 1;
    let flags = record.flags();
    if flags.is_unmapped() {
        stats.unmapped += 1;
        return Ok(None);
    }
    if flags.is_secondary() || flags.is_supplementary() {
        stats.secondary_or_supplementary += 1;
        return Ok(None);
    }
    if flags.is_duplicate() {
        stats.duplicate += 1;
        return Ok(None);
    }
    let mapq = record.mapping_quality().map_or(255, |m| m.get());
    if mapq < p.min_mapq {
        stats.low_mapq += 1;
        return Ok(None);
    }
    let Some(barcode) = get_tag(record, &p.barcode_tag) else {
        stats.missing_barcode += 1;
        return Ok(None);
    };
    if let Some(allow) = &p.barcodes {
        if !allow.contains(&barcode) {
            stats.other_barcode += 1;
            return Ok(None);
        }
    }
    let Some(start) = record.alignment_start().transpose()? else {
        stats.unmapped += 1;
        return Ok(None);
    };
    let start = usize::from(start) as i64 - 1;
    let cigar = cigar_ops(record)?;
    let span = reference_span(&cigar);
    if span <= 0 {
        stats.no_alignment_span += 1;
        return Ok(None);
    }
    stats.kept += 1;
    let hap = haplotype(record, p.hp_tag.as_ref());
    let ps = p
        .ps_tag
        .as_ref()
        .and_then(|t| get_tag(record, t))
        .unwrap_or_else(|| MISSING.to_string());
    let name = record
        .name()
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .unwrap_or_else(|| "*".into());
    let seq = if with_seq {
        record.sequence().iter().collect()
    } else {
        Vec::new()
    };
    Ok(Some((
        (barcode, hap),
        Read {
            name,
            start,
            end: start + span,
            reverse: flags.is_reverse_complemented(),
            mapq,
            ps,
            cigar,
            seq,
        },
    )))
}

struct Parts {
    dir: PathBuf,
}

impl Parts {
    fn path(&self, chrom_idx: usize, kind: &str) -> PathBuf {
        self.dir.join(format!("{chrom_idx:05}.{kind}.tsv"))
    }
    fn create(&self, chrom_idx: usize, kind: &str, on: bool) -> io::Result<Option<BufWriter<File>>> {
        if !on {
            return Ok(None);
        }
        Ok(Some(BufWriter::with_capacity(
            1 << 20,
            File::create(self.path(chrom_idx, kind))?,
        )))
    }
}

type Reader = bam::io::Reader<bgzf::io::MultithreadedReader<File>>;

fn open(bam_path: &Path, workers: usize) -> io::Result<(Reader, noodles_sam::Header)> {
    let file = File::open(bam_path)?;
    let workers = NonZero::new(workers.max(1)).unwrap();
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader.read_header()?;
    Ok((reader, header))
}

fn for_each_record<I, F>(
    bam_path: &Path,
    index: &I,
    chrom: &str,
    workers: usize,
    mut f: F,
) -> io::Result<()>
where
    I: BinningIndex,
    F: FnMut(&bam::Record) -> io::Result<()>,
{
    let (mut reader, header) = open(bam_path, workers)?;
    let region = Region::new(chrom, ..);
    let query = reader.query(&header, index, &region)?;
    for result in query.records() {
        f(&result?)?;
    }
    Ok(())
}

fn process_chrom<I: BinningIndex>(
    bam_path: &Path,
    index: &I,
    chrom: &str,
    chrom_idx: usize,
    chrom_len: i64,
    p: &Params,
    parts: &Parts,
) -> io::Result<ChromOut> {
    let mut out = ChromOut::default();
    let pairs_on = p.write_overlaps || p.write_junctions;

    // optional pre-pass: reads per group, to skip oversized groups
    let mut over: HashSet<GroupKey> = HashSet::new();
    if pairs_on && p.max_reads_per_group > 0 {
        let mut counts: HashMap<GroupKey, u64> = HashMap::new();
        let mut scratch = Stats::default();
        for_each_record(bam_path, index, chrom, p.hts_threads, |rec| {
            if let Some((g, _)) = extract(rec, p, &mut scratch, false)? {
                *counts.entry(g).or_insert(0) += 1;
            }
            Ok(())
        })?;
        for (g, n) in counts {
            if n > p.max_reads_per_group {
                out.skipped.push((g.0.clone(), g.1.clone(), n));
                over.insert(g);
            }
        }
        out.skipped.sort();
    }

    let mut w_ov = parts.create(chrom_idx, "overlaps", p.write_overlaps)?;
    let mut w_jn = parts.create(chrom_idx, "tn5_junctions", p.write_junctions)?;

    let mut active: HashMap<GroupKey, Vec<Read>> = HashMap::new();
    let mut intervals: BTreeMap<GroupKey, Vec<(i64, i64)>> = BTreeMap::new();
    let mut bins: BTreeMap<GroupKey, BTreeMap<i64, u64>> = BTreeMap::new();
    let mut last_start = i64::MIN;
    {
        let mut sink = PairSink::new(
            p,
            chrom,
            w_ov.as_mut().map(|w| w as &mut dyn Write),
            w_jn.as_mut().map(|w| w as &mut dyn Write),
            &mut out.groups,
            &mut out.kmers,
        );
        let stats = &mut out.stats;
        for_each_record(bam_path, index, chrom, p.hts_threads, |rec| {
            let Some((g, read)) = extract(rec, p, stats, pairs_on)? else {
                return Ok(());
            };
            if read.start < last_start {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("{chrom}: BAM is not coordinate-sorted (samtools sort first)"),
                ));
            }
            last_start = read.start;
            sink.groups.entry(g.clone()).or_default().n_reads += 1;
            if p.write_density {
                intervals.entry(g.clone()).or_default().push((read.start, read.end));
                *bins
                    .entry(g.clone())
                    .or_default()
                    .entry(read.start / p.bin_size)
                    .or_insert(0) += 1;
            }
            if pairs_on && !over.contains(&g) {
                let list = active.entry(g.clone()).or_default();
                list.retain(|a| a.end > read.start);
                for a in list.iter() {
                    sink.pair(&g, a, &read)?;
                }
                list.push(read);
            }
            Ok(())
        })?;
    }
    if let Some(w) = w_ov.as_mut() {
        w.flush()?;
    }
    if let Some(w) = w_jn.as_mut() {
        w.flush()?;
    }

    if p.write_density {
        let mut w_seg = parts.create(chrom_idx, "density_segments", true)?.unwrap();
        for (g, iv) in &intervals {
            write_segments(&mut w_seg, g, chrom, iv)?;
        }
        w_seg.flush()?;
        let mut w_bin = parts.create(chrom_idx, "density_bins", true)?.unwrap();
        for (g, counts) in &bins {
            write_bins(&mut w_bin, g, chrom, chrom_len, p.bin_size, counts)?;
        }
        w_bin.flush()?;
    }
    Ok(out)
}

fn gz_create(path: &Path, columns: &[&str]) -> io::Result<GzEncoder<BufWriter<File>>> {
    let mut gz = GzEncoder::new(
        BufWriter::with_capacity(1 << 20, File::create(path)?),
        Compression::new(6),
    );
    writeln!(gz, "{}", columns.join("\t"))?;
    Ok(gz)
}

fn join_parts(parts: &Parts, n: usize, kind: &str, columns: &[&str], dest: &Path) -> io::Result<()> {
    let mut gz = gz_create(dest, columns)?;
    let mut buf = vec![0u8; 1 << 20];
    for i in 0..n {
        let path = parts.path(i, kind);
        if !path.exists() {
            continue;
        }
        let mut f = BufReader::new(File::open(&path)?);
        loop {
            let k = f.read(&mut buf)?;
            if k == 0 {
                break;
            }
            gz.write_all(&buf[..k])?;
        }
        fs::remove_file(&path)?;
    }
    gz.finish()?.flush()
}

/// Run over the given chromosomes (name, length) and write the outputs to
/// `out_dir/{prefix}.{kind}.tsv.gz`.
pub fn run(
    bam_path: &Path,
    out_dir: &Path,
    prefix: &str,
    chroms: &[(String, i64)],
    p: &Params,
    threads: usize,
) -> io::Result<RunOut> {
    fs::create_dir_all(out_dir)?;
    let index = load_index(bam_path)?;
    let parts = Parts {
        dir: out_dir.join(format!(".{prefix}.parts.{}", std::process::id())),
    };
    fs::create_dir_all(&parts.dir)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(io::Error::other)?;
    let results: Vec<io::Result<ChromOut>> = pool.install(|| {
        chroms
            .par_iter()
            .enumerate()
            .map(|(i, (name, len))| match &index {
                Index::Bai(ix) => process_chrom(bam_path, ix, name, i, *len, p, &parts),
                Index::Csi(ix) => process_chrom(bam_path, ix, name, i, *len, p, &parts),
            })
            .collect()
    });

    let mut stats = Stats::default();
    let mut groups: BTreeMap<GroupKey, GroupSummary> = BTreeMap::new();
    let mut kmers: BTreeMap<(String, String, String), u64> = BTreeMap::new();
    let mut skipped = Vec::new();
    let mut first_err = None;
    for ((name, _), r) in chroms.iter().zip(results) {
        match r {
            Ok(o) => {
                stats.add(&o.stats);
                for (g, s) in o.groups {
                    groups.entry(g).or_default().add(&s);
                }
                for (k, n) in o.kmers {
                    *kmers.entry(k).or_insert(0) += n;
                }
                for (b, h, n) in o.skipped {
                    skipped.push((b, h, name.clone(), n));
                }
            }
            Err(e) => {
                first_err.get_or_insert(e);
            }
        }
    }
    if let Some(e) = first_err {
        let _ = fs::remove_dir_all(&parts.dir);
        return Err(e);
    }

    let n = chroms.len();
    let mut files = Vec::new();
    let mut emit = |kind: &str, cols: &[&str]| -> io::Result<()> {
        let dest = out_dir.join(format!("{prefix}.{kind}.tsv.gz"));
        join_parts(&parts, n, kind, cols, &dest)?;
        files.push((kind.to_string(), dest));
        Ok(())
    };
    if p.write_density {
        emit("density_segments", SEGMENT_COLUMNS)?;
        emit("density_bins", BIN_COLUMNS)?;
    }
    if p.write_overlaps {
        emit("overlaps", OVERLAP_COLUMNS)?;
    }
    if p.write_junctions {
        emit("tn5_junctions", OVERLAP_COLUMNS)?;
        let dest = out_dir.join(format!("{prefix}.tn5_summary.tsv.gz"));
        let mut gz = gz_create(&dest, SUMMARY_COLUMNS)?;
        for ((b, h), s) in &groups {
            writeln!(
                gz,
                "{b}\t{h}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.n_reads,
                s.n_overlaps,
                s.n_end_to_start,
                s.n_contained,
                s.n_junctions,
                format_pct(s.n_junctions, s.n_overlaps)
            )?;
        }
        gz.finish()?.flush()?;
        files.push(("tn5_summary".into(), dest));

        let dest = out_dir.join(format!("{prefix}.tn5_kmers.tsv.gz"));
        let mut gz = gz_create(&dest, KMER_COLUMNS)?;
        for ((b, h, k), c) in &kmers {
            writeln!(gz, "{b}\t{h}\t{k}\t{c}")?;
        }
        gz.finish()?.flush()?;
        files.push(("tn5_kmers".into(), dest));
    }
    let _ = fs::remove_dir_all(&parts.dir);
    Ok(RunOut {
        stats,
        skipped,
        files,
    })
}
