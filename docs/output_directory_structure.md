# Output Directory Structure

This page describes the directory structure of `results/` produced by the workflow.

The exact set of directories depends on the target you run and on the configuration in `config/config.yml`.

## Conventions

| Placeholder | Meaning |
| --- | --- |
| `{selected_assembly}` | One of the assemblies listed in `selected_assemblies` in `config/config.yml`: `primary`, `hap1`, or `hap2`. |
| `{assembly_name}` | The value of `assembly_name` in `config/config.yml`. |
| `{organelle}` | One of the organelle types selected by `oatk_organelle`: `mito` or `pltd`. |

- `results/submission/{selected_assembly}/` is created only for assemblies listed in `submission_assemblies`.
- `results/submission/organelle/` depends on `oatk_organelle` and is populated from Oatk outputs.
- `results/ont_reads/` is created only when `ont_reads` is set.
- `results/hic_reads/`, `results/yahs/`, and `results/juicebox/` are created only when both `hic_reads_r1` and `hic_reads_r2` are set.
- Organellar outputs in `results/oatk/` depend on `oatk_organelle` (`mito`, `pltd`, or `mito_and_pltd`).
- Track-specific subdirectories in `results/circos_plot/` depend on `circos_plot_tracks`.

## Downstream Assembly Used by Later Steps

Several downstream steps use a single "final" assembly after contamination removal:

- Without Hi-C reads, the downstream assembly is `results/fcs/assembly/{selected_assembly}/{assembly_name}.fa`.
- With Hi-C reads, the downstream assembly is `results/yahs/assembly/{selected_assembly}/{assembly_name}.fa`.

This affects the inputs for RepeatModeler, RepeatMasker, BRAKER3, submission formatting, and Circos/linear plots.

## `results/braker3/`

BRAKER3 working directories for the RepeatMasker soft-masked downstream assemblies.

- `results/braker3/{selected_assembly}/{assembly_name}/`
  Main tracked outputs include `braker.gff3`, `braker.gtf`, `braker.codingseq`, `braker.aa`, and `augustus_config/`.
  Additional BRAKER3-generated files may also appear in the same directory.

## `results/circos_plot/`

Files used to build the final Circos and linear plots from long contigs of the downstream assembly.

- `results/circos_plot/{selected_assembly}/`
  Final plot PDFs and the long-contig length table.
- `results/circos_plot/gene/{selected_assembly}/`
  BED files and window-based gene coverage tables.
  Created when the `gene` track is enabled.
- `results/circos_plot/repeatmasker_repeat/{selected_assembly}/{repeat_class}/`
  BED files and window-based repeat coverage tables for enabled repeat classes such as `LTR`, `Copia`, `Gypsy`, `LINE`, `SINE`, `DNA_transposon`, and `satellite`.
- `results/circos_plot/tidk_repeat/{selected_assembly}/`
  Window-based telomeric repeat counts from `tidk search`.
  Created when the `tidk` track is enabled.

## `results/downloads/`

Downloaded reference datasets, wrapper scripts, and helper files.

- `busco_downloads/`
  BUSCO lineage datasets used for genome and protein assessments.
- `dfam/`
  Compressed Dfam partitions downloaded before they are unpacked into `results/repeatmasker/dfam/`.
- `fcs/`
  NCBI FCS wrapper scripts, Singularity images, the FCS-GX database (`gxdb/`), and the `.gxdb_checked/` validation directory.
- `oatkdb/`
  Oatk HMM profile database files for the configured lineage and organelle type(s).
- `omamer/`
  The OMAmer database (`LUCA.h5`) used by OMArk.
- `orthodb/`
  OrthoDB protein dataset used by BRAKER3.
- `tidk/`
  The `tidk` database CSV and a cleanup flag file used to restore or remove `~/.local/share/tidk`.

## `results/fcs/`

Assemblies after NCBI FCS processing, plus associated QC outputs.

- `assembly/{selected_assembly}/`
  Final FCS-cleaned assembly copied from `fcs_gx_clean/*.clean.fa`.
- `assembly_long_contigs/{selected_assembly}/`
  Long contigs extracted from the FCS-cleaned assembly.
  Used directly when Hi-C scaffolding is disabled.
- `busco_genome/{selected_assembly}/`
  BUSCO results for the FCS-cleaned assembly.
- `depth/{selected_assembly}/`
  Per-contig depth plots and tables.
  Depth is derived from HiFi read mappings to the Hifiasm assembly and then restricted to contigs present in the FCS-cleaned assembly.
- `fcs_adaptor_screen/{selected_assembly}/`
  FCS-adaptor screening outputs, including reports, logs, JSONL files, YAML arguments, and `cleaned_sequences/`.
- `fcs_adaptor_clean/{selected_assembly}/`
  Assemblies split into cleaned and contaminant FASTA files after adaptor/vector cleaning.
- `fcs_gx_screen/{selected_assembly}/`
  FCS-GX screening reports and taxonomy reports.
- `fcs_gx_clean/{selected_assembly}/`
  Assemblies split into cleaned and contaminant FASTA files after FCS-GX cleaning.
- `gc_content/{selected_assembly}/`
  GC-content tables and plots.
- `length/{selected_assembly}/`
  Contig-length tables and plots.
- `merqury/{selected_assembly}/`
  Merqury QV results.
- `seqkit/{selected_assembly}/`
  SeqKit statistics and exported contig-name lists.
- `tidk/{selected_assembly}/`
  `tidk find`, `tidk explore`, and `tidk search` outputs and plots.

## `results/hic_reads/`

Created only when Hi-C reads are configured.

- `merged/`
  Merged Hi-C FASTQ files for read 1 and read 2.

## `results/hifi_reads/`

HiFi read preprocessing, k-mer analyses, meryl database construction, and read mapping back to the Hifiasm assemblies.

- `raw_reads/`
  FASTQ files converted from the input HiFi BAM files.
- `fastplong/`
  Filtered HiFi reads and per-sample fastplong reports.
- `merged/`
  The merged HiFi FASTQ used in downstream analyses.
- `smudgeplot/`
  FastK outputs and Smudgeplot results.
- `genomescope2/`
  Jellyfish count/histogram files and GenomeScope 2.0 reports and plots.
- `meryl/`
  Meryl databases used by Merqury.
- `map_to_assembly/{selected_assembly}/`
  BAM, BAI, and coverage summary TSV files from mapping merged HiFi reads to Hifiasm assemblies.

## `results/hifiasm/`

Raw Hifiasm outputs, selected assemblies, organelle-screening results, and QC summaries.

- `hifiasm/`
  Raw Hifiasm GFA outputs.
  Without Hi-C input, the primary files are named like `*.asm.bp.*`.
  With Hi-C input, the primary files are named like `*.asm.hic.*`.
- `assembly/{selected_assembly}/`
  FASTA files converted from the selected Hifiasm GFA outputs.
- `assembly_long_contigs/{selected_assembly}/`
  Long contigs extracted from the selected Hifiasm assemblies.
- `busco_genome/{selected_assembly}/`
  BUSCO results for the selected Hifiasm assemblies.
- `depth/{selected_assembly}/`
  Per-contig depth plots and tables.
- `gc_content/{selected_assembly}/`
  GC-content tables and plots.
- `length/{selected_assembly}/`
  Contig-length tables and plots.
- `map_to_organelle/minimap2/{selected_assembly}/`
  PAF alignments of Hifiasm contigs against organelle concatemer references.
- `map_to_organelle/{selected_assembly}/`
  Detected organelle-contig lists, summary tables, and diagnostic plots.
- `map_to_organelle/organelle_contig_depth/{selected_assembly}/`
  Depth plots and tables restricted to the detected organelle contigs.
- `merqury/{selected_assembly}/`
  Merqury QV results.
- `seqkit/{selected_assembly}/`
  SeqKit statistics and exported contig-name lists.
- `tidk/{selected_assembly}/`
  `tidk find`, `tidk explore`, and `tidk search` outputs and plots.

## `results/isoforms/`

Sequence sets copied from BRAKER3 predictions, plus their summary metrics.

- `{selected_assembly}/`
  FASTA files for all predicted CDS, proteins, and transcripts.
- `busco_proteins/{selected_assembly}/`
  BUSCO results for the predicted protein set.
- `seqkit/{selected_assembly}/`
  SeqKit statistics for the predicted protein set.

## `results/juicebox/`

Created only when Hi-C reads are configured.

- `{selected_assembly}/`
  Juicebox-ready files generated from the YaHS scaffolding outputs, including the Hi-C contact map (`*.hic`), assembly annotation file (`*.assembly`), liftover AGP (`*.liftover.agp`), helper AGP (`*.assembly.agp`), link table (`*.txt`), chromosome-size file (`*.chrom.sizes`), and the `juicer pre` log (`*.juicer_pre.log`).
  If the total scaffolded assembly size exceeds 2 Gb, the `*.juicer_pre.log` file also records the Juicebox scale factor that should be set in `Assembly > Set Scale`.

## `results/longest_cds/`

Representative gene models derived from the longest CDS per locus, plus their QC outputs.

- `{selected_assembly}/`
  FASTA files for representative CDS, proteins, and transcripts, plus the filtered representative GFF3 file.
- `busco_proteins/{selected_assembly}/`
  BUSCO results for the representative protein set.
- `omark/{selected_assembly}/`
  OMAmer search output (`*.omamer`) and the OMArk result directory (`*_omark/`).
- `seqkit/{selected_assembly}/`
  SeqKit statistics for the representative protein set.

## `results/oatk/`

Organelle assembly outputs generated from HiFi reads.

- `oatk/`
  Raw Oatk outputs.
  Depending on `oatk_organelle`, this directory can include mitochondrial and/or plastid annotations, GFA files, BED files, and contig FASTA files.
- `concatemer/`
  Concatenated organelle references.
  `*.concatemer.all.fa` is always produced; organelle-specific `*.concatemer.mito.fa` and/or `*.concatemer.pltd.fa` are produced when applicable.
- `seqkit/`
  Organelle-specific SeqKit statistics for mitochondrial and/or plastid contig FASTA files.

## `results/ont_reads/`

Created only when `ont_reads` is set.

- `fastplong/`
  Filtered ONT reads and fastplong reports.

## `results/organelle_removal/`

Assemblies after removing organelle contigs from the selected Hifiasm assemblies, plus associated QC outputs.

- `assembly/{selected_assembly}/`
  Assemblies after organelle-contig removal.
- `assembly_long_contigs/{selected_assembly}/`
  Long contigs extracted from the organelle-removed assemblies.
- `busco_genome/{selected_assembly}/`
  BUSCO results.
- `depth/{selected_assembly}/`
  Per-contig depth plots and tables.
- `gc_content/{selected_assembly}/`
  GC-content tables and plots.
- `length/{selected_assembly}/`
  Contig-length tables and plots.
- `merqury/{selected_assembly}/`
  Merqury QV results.
- `seqkit/{selected_assembly}/`
  SeqKit statistics and exported contig-name lists.
- `tidk/{selected_assembly}/`
  `tidk find`, `tidk explore`, and `tidk search` outputs and plots.

## `results/repeatmasker/`

RepeatMasker inputs and outputs for the downstream assemblies.

- `dfam/`
  Unpacked Dfam FamDB partitions, `dfam_info.txt`, and the lineage-specific repeat FASTA exported from Dfam.
- `library/{selected_assembly}/`
  Merged repeat library combining RepeatModeler and Dfam sequences.
- `{selected_assembly}/`
  Main tracked outputs are the soft-masked FASTA (`*.fa.masked`) and the XML-style annotation table (`*.fa.out.xm`).
  Additional RepeatMasker-generated side files may also appear in this working directory.

## `results/repeatmodeler/`

RepeatModeler database files and de novo repeat-family outputs for the downstream assemblies.

- `{selected_assembly}/`
  BuildDatabase index files (`.nhr`, `.nin`, `.njs`, `.nnd`, `.nni`, `.nog`, `.nsq`, `.translation`) and RepeatModeler outputs such as `*-families.fa`, `*-families.stk`, and `*-rmod.log`.
  Additional RepeatModeler-generated working files such as `RM_*` directories may also appear here.

## `results/rnaseq_reads/`

RNA-seq preprocessing outputs.

- `fastp/`
  Filtered paired-end RNA-seq FASTQ files and the corresponding fastp HTML/JSON reports.

## `results/submission/`

Submission-ready files produced for assemblies listed in `submission_assemblies`, plus organelle submission assets staged from Oatk outputs.

- `{selected_assembly}/`
  Gzipped genome FASTA, isoform CDS/GFF3 files, representative CDS/GFF3 files, and a generated README for submission.
- `organelle/{organelle}/`
  Organelle-specific submission files derived from the corresponding Oatk outputs. Each organelle subdirectory contains a gzipped genome FASTA, a gzipped Oatk annotation text file, and a generated README describing file provenance.

## `results/yahs/`

Created only when Hi-C reads are configured.

- `input/{selected_assembly}/`
  Copies of the FCS-cleaned assemblies together with `samtools faidx` and `bwa-mem2` index files used for scaffolding.
- `alignment/{selected_assembly}/`
  Deduplicated Hi-C BAM files and BAM indices.
- `assembly/{selected_assembly}/`
  Final scaffolded assemblies produced by YaHS.
  These assemblies become the downstream input for RepeatModeler, RepeatMasker, BRAKER3, submission formatting, and plotting when Hi-C is enabled.
- `assembly_long_contigs/{selected_assembly}/`
  Long contigs extracted from the YaHS assemblies for downstream plotting.
- `agp/{selected_assembly}/`
  YaHS AGP files.
- `bin/{selected_assembly}/`
  YaHS binary output files.
- `run/{selected_assembly}/{assembly_name}/`
  YaHS working directory containing raw run outputs copied into the tracked result locations above.
