"""RRecon: read reconstruction for SPLONGGET ONT data.

SPLONGGET fragments genomic DNA with Tn5 mosaic-end tagmentation (as in the
10x Chromium Multiome protocol) while a splint/bridging oligo at the mosaic
ends links short and long fragments from the same original molecule. Because
Tn5 duplicates 9bp of target sequence at each cut, adjacent fragments share a
9bp target-site duplication (TSD) at their junction -- the signal this
package uses to reconnect fragments that were physically split apart.

Modules:
    bam_io      -- generic filtered-read extraction from a BAM file
    overlaps    -- Task 1: per-barcode overlapping-read clustering
"""
