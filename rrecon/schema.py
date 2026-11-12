"""Output tables: file kinds, column orders and column types.

The Rust core (rust/src/core.rs) writes the same columns in the same order;
tests/test_parity.py checks both engines against these lists.

Coordinates are 0-based, half-open ([start, end)), as in BED.
"""

from __future__ import annotations

OVERLAP_COLUMNS = [
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
]
SEGMENT_COLUMNS = ["barcode", "haplotype", "chrom", "start", "end", "depth"]
BIN_COLUMNS = ["barcode", "haplotype", "chrom", "bin_start", "bin_end", "n_reads"]
SUMMARY_COLUMNS = [
    "barcode",
    "haplotype",
    "n_reads",
    "n_overlaps",
    "n_end_to_start",
    "n_contained",
    "n_junctions",
    "pct_junctions",
]
KMER_COLUMNS = ["barcode", "haplotype", "junction_seq", "count"]

COLUMNS = {
    "density_segments": SEGMENT_COLUMNS,
    "density_bins": BIN_COLUMNS,
    "overlaps": OVERLAP_COLUMNS,
    "tn5_junctions": OVERLAP_COLUMNS,
    "tn5_summary": SUMMARY_COLUMNS,
    "tn5_kmers": KMER_COLUMNS,
}

# which files each output group produces
OUTPUT_KINDS = {
    "density": ["density_segments", "density_bins"],
    "overlaps": ["overlaps"],
    "junctions": ["tn5_junctions", "tn5_summary", "tn5_kmers"],
}

INT_COLUMNS = {
    "read1_start",
    "read1_end",
    "read1_mapq",
    "read2_start",
    "read2_end",
    "read2_mapq",
    "overlap_start",
    "overlap_end",
    "overlap_bp",
    "n_mismatches",
    "start",
    "end",
    "depth",
    "bin_start",
    "bin_end",
    "n_reads",
    "n_overlaps",
    "n_end_to_start",
    "n_contained",
    "n_junctions",
    "count",
}
BOOL_COLUMNS = {"has_indel_in_overlap", "seq_match"}
FLOAT_COLUMNS = {"pct_junctions"}

HAPLOTYPES = ("1", "2", "untagged")
UNTAGGED = "untagged"
MISSING = "."
