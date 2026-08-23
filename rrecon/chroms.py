"""Canonical chromosome definitions shared across rrecon modules.

"Canonical" here means the standard human autosomes (1-22) and sex
chromosomes (X, Y) -- explicitly excluding chrM/MT, alt haplotypes, decoy
sequences, and unplaced/unlocalized scaffolds that show up in reference
genomes like GRCh38's full analysis set. Restricting to these upfront
matters for two reasons: correctness (alt/decoy contigs attract ambiguous,
often multi-mapping reads that would pollute overlap counts) and speed
(skipping them at the chromosome-dispatch level in find_overlaps.py means
they're never scanned at all, not just filtered out after the fact).
"""

_AUTOSOMES = [str(i) for i in range(1, 23)]
_SEX_CHROMS = ["X", "Y"]
_BARE_NAMES = _AUTOSOMES + _SEX_CHROMS

# Covers both naming conventions a reference might use ("chr1" vs "1") --
# we haven't confirmed which your BAM follows, so both are accepted.
DEFAULT_CANONICAL_CHROMS = frozenset(_BARE_NAMES) | frozenset(f"chr{n}" for n in _BARE_NAMES)


def is_canonical_chrom(chrom: str, canonical=DEFAULT_CANONICAL_CHROMS) -> bool:
    """True if `chrom` is one of the standard autosomes/sex chromosomes."""
    return chrom in canonical
