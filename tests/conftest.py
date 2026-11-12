import pytest

from rrecon.synthetic import case_bam


@pytest.fixture(scope="session")
def bam(tmp_path_factory):
    """Synthetic BAM with every edge case (rrecon.synthetic.case_reads)."""
    return case_bam(tmp_path_factory.mktemp("bam") / "cases.bam")
