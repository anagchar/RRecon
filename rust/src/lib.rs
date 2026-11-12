//! rrecon core: read density, pairwise read overlaps and Tn5 junctions from a
//! coordinate-sorted, indexed BAM. Exposed to Python as `rrecon._core`.

pub mod align;
pub mod core;
pub mod scan;

#[cfg(feature = "python")]
mod python {
    use std::collections::HashSet;
    use std::path::Path;

    use pyo3::exceptions::{PyIOError, PyValueError};
    use pyo3::prelude::*;

    use crate::core::Params;

    fn tag(s: &str) -> PyResult<[u8; 2]> {
        let b = s.as_bytes();
        if b.len() != 2 {
            return Err(PyValueError::new_err(format!("BAM tag must be 2 characters, got {s:?}")));
        }
        Ok([b[0], b[1]])
    }

    type RunResult = (
        Vec<(String, u64)>,
        Vec<(String, String, String, u64)>,
        Vec<(String, String)>,
    );

    /// Run the scan; returns (stats, skipped_groups, files).
    #[pyfunction]
    #[pyo3(signature = (
        bam, out_dir, prefix, chroms, *, barcode_tag, hp_tag, ps_tag, min_mapq, barcodes,
        bin_size, junction_size, overlap_tolerance, max_mismatches, max_reads_per_group,
        max_seq_bp, write_density, write_overlaps, write_junctions, threads, hts_threads
    ))]
    #[allow(clippy::too_many_arguments)]
    fn run(
        py: Python<'_>,
        bam: String,
        out_dir: String,
        prefix: String,
        chroms: Vec<(String, i64)>,
        barcode_tag: String,
        hp_tag: Option<String>,
        ps_tag: Option<String>,
        min_mapq: u8,
        barcodes: Option<Vec<String>>,
        bin_size: i64,
        junction_size: i64,
        overlap_tolerance: i64,
        max_mismatches: i64,
        max_reads_per_group: u64,
        max_seq_bp: i64,
        write_density: bool,
        write_overlaps: bool,
        write_junctions: bool,
        threads: usize,
        hts_threads: usize,
    ) -> PyResult<RunResult> {
        if bin_size <= 0 {
            return Err(PyValueError::new_err("bin_size must be positive"));
        }
        let p = Params {
            barcode_tag: tag(&barcode_tag)?,
            hp_tag: hp_tag.as_deref().map(tag).transpose()?,
            ps_tag: ps_tag.as_deref().map(tag).transpose()?,
            min_mapq,
            barcodes: barcodes.map(|v| v.into_iter().collect::<HashSet<_>>()),
            bin_size,
            junction_size,
            overlap_tolerance,
            max_mismatches,
            max_reads_per_group,
            max_seq_bp,
            write_density,
            write_overlaps,
            write_junctions,
            hts_threads,
        };
        let out = py
            .detach(|| {
                crate::scan::run(
                    Path::new(&bam),
                    Path::new(&out_dir),
                    &prefix,
                    &chroms,
                    &p,
                    threads,
                )
            })
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok((
            out.stats.to_vec(),
            out.skipped,
            out.files
                .into_iter()
                .map(|(k, v)| (k, v.to_string_lossy().into_owned()))
                .collect(),
        ))
    }

    #[pymodule]
    fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(run, m)?)?;
        m.add("__version__", env!("CARGO_PKG_VERSION"))?;
        Ok(())
    }
}
