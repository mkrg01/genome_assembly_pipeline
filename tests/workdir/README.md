This directory contains tiny placeholder inputs used only for Snakemake dry-run checks in CI.

The files under `raw_data/` are not valid sequencing datasets and must not be used for real workflow execution.

The VARUS annotation check in CI copies only the configuration and assembly FASTA
to a temporary work directory, verifying that nuclear annotation in VARUS mode
does not depend on local RNA-seq. Dry-runs do not query SRA or download reads.
