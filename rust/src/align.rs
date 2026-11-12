//! CIGAR-aware extraction of the bases a read places on a reference window.
//!
//! The Python fallback (`rrecon/_pyengine.py::aligned_window`) implements the
//! same rules; the parity tests check that both give identical output.

/// CIGAR operation kinds, in BAM code order (M I D N S H P = X).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpKind {
    Match,
    Ins,
    Del,
    Skip,
    Soft,
    Hard,
    Pad,
    SeqMatch,
    SeqMismatch,
}

impl OpKind {
    pub fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            0 => Self::Match,
            1 => Self::Ins,
            2 => Self::Del,
            3 => Self::Skip,
            4 => Self::Soft,
            5 => Self::Hard,
            6 => Self::Pad,
            7 => Self::SeqMatch,
            8 => Self::SeqMismatch,
            _ => return None,
        })
    }

    pub fn consumes_reference(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Del | Self::Skip | Self::SeqMatch | Self::SeqMismatch
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CigarOp {
    pub kind: OpKind,
    pub len: u32,
}

impl CigarOp {
    pub fn new(kind: OpKind, len: u32) -> Self {
        Self { kind, len }
    }
}

/// Number of reference bases the alignment spans.
pub fn reference_span(cigar: &[CigarOp]) -> i64 {
    cigar
        .iter()
        .filter(|op| op.kind.consumes_reference())
        .map(|op| op.len as i64)
        .sum()
}

/// The read's bases on reference window `[ws, we)`, one per reference position.
///
/// - aligned bases (M, =, X) are copied, upper-cased; a read without a stored
///   sequence gives `N`
/// - reference positions the read deletes or skips (D, N), or does not cover,
///   are `-`
/// - soft clips are skipped, hard clips and padding ignored
///
/// Also returns whether the read has an insertion strictly inside the window,
/// i.e. between two reference positions that both lie in `[ws, we)`.
pub fn aligned_window(
    ref_start: i64,
    cigar: &[CigarOp],
    seq: &[u8],
    ws: i64,
    we: i64,
) -> (Vec<u8>, bool) {
    let n = (we - ws).max(0) as usize;
    let mut out = vec![b'-'; n];
    let mut insertion = false;
    let mut rpos = ref_start;
    let mut qpos: usize = 0;
    for op in cigar {
        if rpos >= we {
            break;
        }
        let len = op.len as i64;
        match op.kind {
            OpKind::Match | OpKind::SeqMatch | OpKind::SeqMismatch => {
                let lo = rpos.max(ws);
                let hi = (rpos + len).min(we);
                for p in lo..hi {
                    let q = qpos + (p - rpos) as usize;
                    out[(p - ws) as usize] = base_at(seq, q);
                }
                rpos += len;
                qpos += op.len as usize;
            }
            OpKind::Ins => {
                if rpos > ws && rpos < we {
                    insertion = true;
                }
                qpos += op.len as usize;
            }
            OpKind::Del | OpKind::Skip => rpos += len,
            OpKind::Soft => qpos += op.len as usize,
            OpKind::Hard | OpKind::Pad => {}
        }
    }
    (out, insertion)
}

fn base_at(seq: &[u8], q: usize) -> u8 {
    seq.get(q).map_or(b'N', |b| b.to_ascii_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;
    use OpKind::*;

    fn s(v: &[u8]) -> String {
        String::from_utf8(v.to_vec()).unwrap()
    }

    #[test]
    fn plain_match() {
        let c = [CigarOp::new(Match, 10)];
        let (w, ins) = aligned_window(100, &c, b"ACGTACGTAC", 103, 107);
        assert_eq!(s(&w), "TACG");
        assert!(!ins);
    }

    #[test]
    fn soft_clip_is_skipped() {
        let c = [CigarOp::new(Soft, 3), CigarOp::new(Match, 5)];
        let (w, _) = aligned_window(10, &c, b"nnnACGTA", 10, 15);
        assert_eq!(s(&w), "ACGTA");
    }

    #[test]
    fn deletion_inside_window() {
        // ref 0..6, read ACG - - T  (2D)
        let c = [CigarOp::new(Match, 3), CigarOp::new(Del, 2), CigarOp::new(Match, 1)];
        let (w, ins) = aligned_window(0, &c, b"ACGT", 0, 6);
        assert_eq!(s(&w), "ACG--T");
        assert!(!ins);
    }

    #[test]
    fn insertion_inside_and_at_edges() {
        let c = [CigarOp::new(Match, 3), CigarOp::new(Ins, 2), CigarOp::new(Match, 3)];
        let seq = b"ACGTTACG";
        // insertion sits between ref 2 and 3
        assert!(aligned_window(0, &c, seq, 0, 6).1);
        let (w, _) = aligned_window(0, &c, seq, 0, 6);
        assert_eq!(s(&w), "ACGACG");
        // window starting at 3: the insertion is on its left edge, not inside
        assert!(!aligned_window(0, &c, seq, 3, 6).1);
        // window ending at 3: on its right edge
        assert!(!aligned_window(0, &c, seq, 0, 3).1);
    }

    #[test]
    fn missing_sequence_gives_n() {
        let c = [CigarOp::new(Match, 4)];
        let (w, _) = aligned_window(0, &c, b"", 1, 3);
        assert_eq!(s(&w), "NN");
    }

    #[test]
    fn span_counts_reference_ops() {
        let c = [
            CigarOp::new(Soft, 5),
            CigarOp::new(Match, 10),
            CigarOp::new(Ins, 3),
            CigarOp::new(Del, 2),
            CigarOp::new(SeqMatch, 4),
        ];
        assert_eq!(reference_span(&c), 16);
    }
}
