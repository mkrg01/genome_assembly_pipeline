# BRAKER4 execution and migration

Gene prediction uses the official BRAKER4 **v0.5.0-beta** release, commit
`c72618baa4631db4a4128391a40c8ba7392df17c`. The source archive is verified with
SHA256 `c4f22285595c1945fb97a0a4ef25a047e12e097a1eb595c09b230771c7a3fe42`.
The release's main tool image is `teambraker/braker3:v3.0.10`; the image retains
its upstream name, but this workflow never invokes the old `braker.pl` controller.
AGAT is pinned to `1.4.1--pl5321hdfd78af_0`.

## Running annotation

Run from the repository root, using the existing input/configuration conventions:

```bash
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" \
  --cores 48 --resources mem_mb=120000 gene_prediction_all
```

For external assemblies, add `-s workflow/Snakefile.annotation` and set
`external_assembly` as described in the main README. Per-assembly resources are
defined directly in `rule braker4` in `workflow/rules/gene_prediction.smk`:
`threads: 48` and `resources: mem_mb=120000`; adjust these rule values to the
resources available. Snakemake caps the child CPU allocation to the parent job's
allocation. Memory is a scheduling budget, not an operating-system memory limit.
Use a total parent `--resources mem_mb=...` budget when several assemblies can run
concurrently. The existing SLURM and SHIROKANE job scripts can launch this workflow
within their allocated node; the child uses the local executor.

The BRAKER4 controller explicitly opts out of the global pipeline container. Its
host environment, `workflow/envs/braker4.yml`, pins Snakemake 9.8.0 and pandas 2.2.3.
It requires **Conda and Apptainer on the host**, and both deployment methods must
be enabled. This host-only environment does not require rebuilding the existing
pipeline Docker image. The controller launches the tool containers itself, with
explicit bindings for its inputs, source, database, cache and temporary directory.
Parent cluster profiles and inherited `BRAKER4_*` overrides are not used by the
child; generated configuration determines its behavior.

The source, controller environment, images, OrthoDB and BUSCO lineage must be
downloaded before an offline compute-node run. `--conda-create-envs-only` on the
parent can prepare the host environment; the child images are downloaded by the
child, not by the parent's `--singularity-pull` operation. Its shared image cache
is `results/downloads/braker4/containers/`. Offline execution applies to local
RNA-seq mode; VARUS requires network access during sampling.

BRAKER4's upstream shell commands require input and working paths without spaces
or shell metacharacters. The wrapper checks this before prediction. Its sample
name, `{organism_name}_{selected_assembly}`, must start with a letter and contain
only letters, digits and underscores.

## RNA-seq source selection

The default is `braker4_rnaseq_source: "local"`: all detected local RNA-seq pairs
are processed with fastp and supplied to BRAKER4. To enable BRAKER4's automated
RNA-seq sampling, set the following in `config/config.yml`:

```yaml
braker4_rnaseq_source: "varus"
```

The genus and species come directly from `organism_name` (for example,
`Arabidopsis_thaliana`). VARUS requires this to be a binomial scientific name,
without cultivar, accession or subspecies components.
The query populates BRAKER4's `varus_genus` and `varus_species` sample fields and
enables `use_varus`. VARUS uses the pinned upstream image
`katharinahoff/varus-notebook:v0.0.6` to sample SRA reads and produce the RNA-seq BAM
for ETP annotation with the existing OrthoDB proteins.

VARUS mode does not require local RNA-seq for `gene_prediction_all`,
`circos_plot_all`, or the external-assembly annotation workflow. Local libraries
are not combined with VARUS evidence for nuclear prediction. The selected mode
is explicit; missing or unsuitable SRA evidence does not trigger a fallback to
local reads or protein-only prediction. The compute node needs internet access
to SRA during sampling, and suitable public RNA-seq must exist for the query.
Download volume and runtime depend on the available data and upstream sampling
settings; this integration does not expose a download-size limit.

The full assembly workflow's organelle RNA-editing curation still uses local
paired-end RNA-seq mapped against the combined nuclear/organelle reference.
BRAKER4's nuclear VARUS BAM cannot replace this evidence. Supply local reads for
that stage, choose a nuclear annotation target, or explicitly disable organelle
annotation with `organelle_annotation: null` if it is not needed.

Sampling runs independently for each selected assembly. Preserve the work
directory and its native `output/{sample}/varus/` data to resume without sampling
again. The native results directory also retains `varus_runlist.tsv` and
`varus_stats.txt`, recording accessions and download statistics. `run.json`
records the selected RNA-seq source, genus/species query, and native results path.
A fresh run may select different evidence as the public SRA collection changes.

## Evidence and annotation settings

- The RepeatMasker soft-masked genome populates both `genome` and `genome_masked`.
  BRAKER4's repeat masking is therefore skipped, preserving the existing mask.
- In local mode, all fastp-filtered RNA-seq pairs populate `fastq_r1` and `fastq_r2` in matching
  order. Unique local aliases prevent library names containing dots from colliding
  under upstream filename parsing. Local sample names are never treated as SRA IDs.
- The existing OrthoDB proteins populate `protein_fasta`, selecting ETP mode.
- `busco_lineage_dataset` must use **odb12**. The parent downloads the lineage;
  compleasm reuses `results/downloads/busco_downloads/lineages/`.
- `workflow/config/braker4.ini` retains AUGUSTUS optimization and compleasm rescue,
  excludes compleasm-derived hints from the prediction evidence, and disables
  duplicate BUSCO/OMArk assessments inside BRAKER4. The existing downstream QC
  remains in place. ncRNA and FANTASIA are disabled.

In this release the GFF3 is generated from the UTR-decorated GTF when transcript
evidence is available. Thus transcripts extracted by the parent from this GFF3
can include supported UTRs. `braker.gtf` remains the coding-model GTF. Both files
and the CDS/protein FASTAs must contain matching transcript IDs before export.
This integration supports local paired-end or VARUS short-read ETP annotation; adding Iso-Seq
inputs would require a separate extension of the input configuration.

## Outputs and restart

The parent consumes validated, decompressed files in
`results/braker4/{selected_assembly}/{organism_name}/`:
`braker.gff3`, `braker.gtf`, `braker.codingseq`, `braker.aa`, and `run.json`.
`run.json` records the source version, input identities, settings hash, child
command, successful work directory, and gene/transcript counts.

The child runs in `work/{run_id}/` below that directory. Its ID is derived from
the source manifest, annotation template, wrapper, sample, lineage, RNA-seq mode,
VARUS query (when used), and input
paths/sizes/modification times. Changing inputs or annotation settings starts a
new directory; changing only CPU/memory allocation allows a retry to resume.
Input content changed in place while preserving both size and modification time
is not detected by this identity scheme.

The child runs with `--rerun-incomplete`, and all intermediates and AUGUSTUS
training files are retained. Rerunning the parent after a failure reuses completed
child steps. The parent never declares the child working directory as a
`directory()` output, which would allow Snakemake to remove it before a retry.
Inspect `logs/braker4_*.out` and `logs/braker4_*.err` for the failed run's path.
Work directories can be large; removing one discards its restart state.

The exporter checks for completed collection, nonempty/unique sequence IDs,
matching CDS/protein/GTF/GFF3 transcript sets, and valid gene/mRNA/CDS relationships.
It stages all four files before publishing them. Longest-CDS selection uses GFF3
parent genes, and submission ID conversion also supports independent transcript
IDs introduced by annotation formatting.

## Migrating an existing analysis

The BRAKER3 execution branch has been removed. Existing `results/braker3/` files
are left on disk but are no longer workflow inputs. New runs use `results/braker4/`.
Archive the old QC and release files first if a comparison is needed: the downstream
`isoforms/`, `longest_cds/`, `release/`, and plot paths are still shared with the
previous pipeline. For a new published annotation, use a new `genome_version`.

Inspect the normal `gene_prediction_all` dry-run before executing; upstream
assemblies, masking and RNA-seq preprocessing can be reused. For the first
migration, use Snakemake's normal rerun triggers rather than restricting them to
mtime, so changed input paths and code are detected.

Before releasing annotations, compare a representative full-genome run against
the archived gene counts, BUSCO/OMArk results and RNA-seq support. Tiny CI fixtures
verify workflow wiring, not GeneMark training or biological annotation quality.

## Validation for developers

```bash
python -m unittest discover -s tests -p 'test_*.py'
snakemake --lint
snakemake --directory tests/workdir --cores 2 --dry-run all
```

To additionally dry-run the real pinned BRAKER4 DAG, use the controller environment
and set `BRAKER4_TEST_SOURCE` to its extracted source directory before running
`python -m unittest discover -s tests -p 'test_braker4.py'`. These optional tests
cover local ETP with two synthetic RNA-seq libraries and VARUS ETP with no local
reads, using a pre-masked genome, proteins and an odb12 lineage placeholder.
Apptainer must be installed; tool images and SRA sampling are not executed.

References: [release](https://github.com/Gaius-Augustus/BRAKER4/releases/tag/v0.5.0-beta),
[release configuration](https://github.com/Gaius-Augustus/BRAKER4/blob/v0.5.0-beta/config.ini.example),
[release workflow](https://github.com/Gaius-Augustus/BRAKER4/blob/v0.5.0-beta/Snakefile),
[VARUS rule](https://github.com/Gaius-Augustus/BRAKER4/blob/v0.5.0-beta/rules/preprocessing/run_varus.smk).
