import argparse
import configparser
import csv
import gzip
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflow/scripts"))
import braker4_compleasm
import download_braker4
import filter_gff_by_fasta_ids
import rename_submission_gene_models
import run_braker4
import select_longest_cds


GFF3 = """##gff-version 3
chr1\tBRAKER4\tgene\t1\t100\t.\t+\t.\tID=geneA;
chr1\tBRAKER4\tmRNA\t1\t100\t.\t+\t.\tID=transcriptA;Parent=geneA;
chr1\tBRAKER4\tCDS\t1\t6\t.\t+\t0\tID=cdsA;Parent=transcriptA;
chr1\tBRAKER4\tmRNA\t1\t100\t.\t+\t.\tID=transcriptB;Parent=geneA;
chr1\tBRAKER4\tCDS\t1\t9\t.\t+\t0\tID=cdsB;Parent=transcriptB;
chr1\tBRAKER4\tgene\t200\t300\t.\t+\t.\tID=geneB;
chr1\tBRAKER4\tmRNA\t200\t300\t.\t+\t.\tID=independent;Parent=geneB;
chr1\tBRAKER4\tCDS\t200\t208\t.\t+\t0\tID=cdsC;Parent=independent;
"""
CDS = ">transcriptA\nATGAAA\n>transcriptB\nATGAAAAAA\n>independent\nATGAAAAAA\n"
AA = ">transcriptA\nMK\n>transcriptB\nMKK\n>independent\nMKK\n"
GTF = """chr1\tBRAKER4\tCDS\t1\t6\t.\t+\t0\tgene_id "geneA"; transcript_id "transcriptA";
chr1\tBRAKER4\tCDS\t1\t9\t.\t+\t0\tgene_id "geneA"; transcript_id "transcriptB";
chr1\tBRAKER4\tCDS\t200\t208\t.\t+\t0\tgene_id "geneB"; transcript_id "independent";
"""
NATIVE_GFF3 = GFF3.replace("\tmRNA\t", "\ttranscript\t") + (
    "chr1\tAGAT\tmRNA\t1\t100\t.\t+\t.\tID=agat-mrna-1;Parent=geneA;transcript_id=transcriptB;\n"
    "chr1\tStringTie\tthree_prime_UTR\t10\t12\t.\t+\t.\tParent=transcriptB;\n"
)


class Braker4Test(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        source = self.root / "source"
        source.mkdir()
        (source / "Snakefile").write_text("rule all:\n    input: []\n")
        manifest = self.root / "source.json"
        manifest.write_text(json.dumps({"version": "v0.8.1-beta", "sha256": "test"}))
        for name in ("genome.fa", "proteins.fa"):
            (self.root / name).write_text(">chr1\nATGAAA\n")
        for lib in ("tissue.rep1", "tissue.rep2"):
            for mate in (1, 2):
                (self.root / f"{lib}_{mate}.fastq").write_text("@r\nATG\n+\nIII\n")
        lineage = self.root / "busco/lineages/embryophyta_odb12"
        lineage.mkdir(parents=True)
        (lineage / "dataset.cfg").write_text("name=embryophyta_odb12\n")
        self.args = argparse.Namespace(
            source=source, source_manifest=manifest,
            template=REPO / "workflow/config/braker4.ini",
            sample="Testus_example_primary", genome=self.root / "genome.fa",
            proteins=self.root / "proteins.fa", lineage="embryophyta_odb12",
            rnaseq_source="local", varus_genus="", varus_species="",
            rnaseq_r1=[self.root / f"tissue.rep{i}_1.fastq" for i in (1, 2)],
            rnaseq_r2=[self.root / f"tissue.rep{i}_2.fastq" for i in (1, 2)],
            busco_downloads=self.root / "busco", threads=2, mem_mb=32000,
            container_cache=self.root / "containers", output_dir=self.root / "results",
            dry_run=False,
        )

    def outputs(self, run):
        results = run / "output" / self.args.sample / "results"
        results.mkdir(parents=True)
        for name, contents in zip(run_braker4.CORE_FILES, (GTF, GFF3, CDS, AA)):
            with gzip.open(results / (name + ".gz"), "wt") as handle:
                handle.write(contents)
        (results / ".done").touch()
        sample_dir = results.parent
        (sample_dir / "best_by_compleasm.log").write_text("STATUS: OK: original gene set kept\n")
        summary = "S:90.0%\nD:1.0%\nF:2.0%\nM:7.0%\n"
        for name in ("genemark", "augustus", "braker"):
            directory = sample_dir / "best_by_compleasm_tmp/single" / name
            directory.mkdir(parents=True)
            (directory / "summary.txt").write_text(summary)
        directory = sample_dir / "compleasm_proteins"
        directory.mkdir()
        (directory / "summary.txt").write_text(summary)
        return results

    def test_etp_configuration_and_unique_library_aliases(self):
        with mock.patch.dict(os.environ, {"BRAKER4_USE_COMPLEASM_HINTS": "1"}):
            run, _, cmd, env = run_braker4.prepare_run(self.args)
        with (run / "samples.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["genome"], row["genome_masked"])
        self.assertEqual(row["protein_fasta"], str(self.args.proteins))
        self.assertEqual("", row["sra_ids"])
        self.assertEqual("", row["bam_files"])
        r1 = row["fastq_r1"].split(":")
        r2 = row["fastq_r2"].split(":")
        self.assertEqual(2, len(r1))
        self.assertEqual(2, len({Path(p).name.split('.')[0] for p in r1}))
        for actual, expected in zip(r1 + r2, self.args.rnaseq_r1 + self.args.rnaseq_r2):
            self.assertEqual(Path(actual).resolve(), expected)
        ini = configparser.ConfigParser()
        ini.read(run / "config.ini")
        self.assertEqual(str(run / "compleasm_downloads"), ini["paths"]["compleasm_download_path"])
        self.assertEqual(2, ini.getint("SLURM_ARGS", "cpus_per_task"))
        self.assertEqual(32000, ini.getint("SLURM_ARGS", "mem_of_node"))
        self.assertTrue(ini.getboolean("PARAMS", "skip_busco"))
        self.assertTrue(ini.getboolean("PARAMS", "no_cleanup"))
        self.assertFalse(ini.getboolean("PARAMS", "run_omark"))
        self.assertFalse(ini.getboolean("PARAMS", "use_varus"))
        self.assertEqual("", row["varus_genus"])
        self.assertEqual("", row["varus_species"])
        self.assertNotIn("BRAKER4_USE_COMPLEASM_HINTS", env)
        self.assertEqual("local", cmd[cmd.index("--executor") + 1])
        self.assertIn("mem_mb=32000", cmd)
        self.assertEqual("none", cmd[cmd.index("--workflow-profile") + 1])
        self.assertIn(str(self.root), cmd[cmd.index("--singularity-args") + 1])

    def use_varus(self):
        self.args.rnaseq_source = "varus"
        self.args.varus_genus = "Arabidopsis"
        self.args.varus_species = "thaliana"
        self.args.rnaseq_r1 = []
        self.args.rnaseq_r2 = []

    def test_varus_without_local_reads_and_query_changes(self):
        local, _, _, _ = run_braker4.prepare_run(self.args)
        self.use_varus()
        for path in self.root.glob("*.fastq"):
            path.unlink()
        run, identity, _, _ = run_braker4.prepare_run(self.args)
        self.assertNotEqual(local, run)
        with (run / "samples.csv").open() as handle:
            row = next(csv.DictReader(handle))
        for field in ("fastq_r1", "fastq_r2", "sra_ids", "bam_files"):
            self.assertEqual("", row[field])
        self.assertEqual("Arabidopsis", row["varus_genus"])
        self.assertEqual("thaliana", row["varus_species"])
        self.assertEqual("varus", identity["rnaseq_source"])
        self.assertEqual([], identity["rnaseq_r1"])
        ini = configparser.ConfigParser()
        ini.read(run / "config.ini")
        self.assertTrue(ini.getboolean("PARAMS", "use_varus"))
        self.assertEqual("docker://katharinahoff/varus-notebook:v0.0.6", ini["containers"]["varus_image"])
        self.assertEqual(run, run_braker4.prepare_run(self.args)[0])
        self.args.varus_species = "lyrata"
        self.assertNotEqual(run, run_braker4.prepare_run(self.args)[0])

    def test_invalid_or_mixed_rnaseq_modes_are_rejected(self):
        self.args.varus_genus = "Arabidopsis"
        with self.assertRaisesRegex(ValueError, "cannot be supplied"):
            run_braker4.prepare_run(self.args)
        self.use_varus()
        self.args.rnaseq_r1 = [self.root / "tissue.rep1_1.fastq"]
        with self.assertRaisesRegex(ValueError, "do not also supply"):
            run_braker4.prepare_run(self.args)
        self.args.rnaseq_r1 = []
        for invalid in ("", "thaliana;echo", "thaliana extra"):
            self.args.varus_species = invalid
            with self.assertRaisesRegex(ValueError, "valid genus and species"):
                run_braker4.prepare_run(self.args)
        self.args.rnaseq_source = "local"
        with self.assertRaisesRegex(ValueError, "matching, non-empty"):
            run_braker4.prepare_run(self.args)
        self.args.rnaseq_source = "automatic"
        with self.assertRaisesRegex(ValueError, "must be 'local' or 'varus'"):
            run_braker4.prepare_run(self.args)

    def test_retries_resume_but_changed_evidence_gets_a_new_workdir(self):
        first, _, _, _ = run_braker4.prepare_run(self.args)
        checkpoint = first / "completed-training"
        checkpoint.touch()
        self.args.threads = 1
        second, _, _, _ = run_braker4.prepare_run(self.args)
        self.assertEqual(first, second)
        self.assertTrue(checkpoint.exists())
        self.args.rnaseq_r1.pop()
        self.args.rnaseq_r2.pop()
        third, _, _, _ = run_braker4.prepare_run(self.args)
        self.assertNotEqual(first, third)
        self.assertTrue(checkpoint.exists())

    def test_assemblies_have_separate_training_directories(self):
        first, _, _, _ = run_braker4.prepare_run(self.args)
        self.args.sample = "Testus_example_hap1"
        second, _, _, _ = run_braker4.prepare_run(self.args)
        self.assertNotEqual(first, second)

    def test_missing_pair_or_lineage_is_rejected(self):
        self.args.rnaseq_r2.pop()
        with self.assertRaisesRegex(ValueError, "matching"):
            run_braker4.prepare_run(self.args)
        self.args.rnaseq_r1.pop()
        self.args.lineage = "embryophyta_odb10"
        with self.assertRaisesRegex(ValueError, "odb12"):
            run_braker4.prepare_run(self.args)

    def test_exports_compressed_results_through_submission_processing(self):
        run, identity, cmd, _ = run_braker4.prepare_run(self.args)
        results = self.outputs(run)
        with gzip.open(results / "braker.gff3.gz", "wt") as handle:
            handle.write(NATIVE_GFF3)
        run_braker4.export_results(run, self.args.output_dir, self.args.sample, identity, cmd)
        exported = self.args.output_dir
        self.assertEqual({"genes": 2, "transcripts": 3}, json.loads((exported / "run.json").read_text())["counts"])
        gff3 = (exported / "braker.gff3").read_text()
        self.assertNotIn("agat-mrna-1", gff3)
        self.assertIn("\tthree_prime_UTR\t10\t12\t.\t+\t.\tParent=transcriptB;", gff3)
        representative = exported / "representative.fa"
        select_longest_cds.select_longest(exported / "braker.codingseq", exported / "braker.gff3", representative)
        self.assertEqual({"transcriptB", "independent"}, run_braker4.fasta_ids(representative))
        filtered = exported / "representative.gff3"
        with mock.patch.object(sys, "argv", ["filter", "--fasta", str(representative),
                                             "--gff3", str(exported / "braker.gff3"), "--output", str(filtered)]):
            filter_gff_by_fasta_ids.main()
        genes, transcripts = rename_submission_gene_models.build_identifier_maps(exported / "braker.gff3", "Tesex")
        rename_submission_gene_models.rewrite_gff3(filtered, exported / "submission.gff3", genes, transcripts)
        rename_submission_gene_models.rewrite_fasta(representative, exported / "submission.fa", transcripts)
        self.assertEqual({"Tesex_000001.t2", "Tesex_000002.t1"}, run_braker4.fasta_ids(exported / "submission.fa"))
        self.assertIn("Parent=Tesex_000001.t2", (exported / "submission.gff3").read_text())

    def test_referenced_or_inconsistent_agat_aliases_are_rejected(self):
        run, identity, cmd, _ = run_braker4.prepare_run(self.args)
        results = self.outputs(run)
        for gff3 in (
            NATIVE_GFF3.replace("Parent=transcriptB;\n", "Parent=agat-mrna-1;\n"),
            NATIVE_GFF3.replace("Parent=geneA;transcript_id=transcriptB", "Parent=geneB;transcript_id=transcriptB"),
        ):
            with self.subTest(gff3=gff3.splitlines()[-2:]):
                with gzip.open(results / "braker.gff3.gz", "wt") as handle:
                    handle.write(gff3)
                with self.assertRaisesRegex(ValueError, "AGAT mRNA alias"):
                    run_braker4.export_results(run, self.args.output_dir, self.args.sample, identity, cmd)
                self.assertFalse((self.args.output_dir / "run.json").exists())

    def test_reexports_completed_run_with_same_inputs_from_an_older_runner(self):
        run, identity, _, _ = run_braker4.prepare_run(self.args)
        previous = {**identity, "runner_sha256": "a" * 64}
        older_run = run.with_name(run_braker4.run_name(previous))
        older_run.mkdir()
        (older_run / "inputs.json").write_text(json.dumps(previous))
        self.outputs(older_run)
        self.args.export_completed_run = older_run

        run_braker4.export_completed_run(self.args)
        manifest = json.loads((self.args.output_dir / "run.json").read_text())
        self.assertEqual(previous["runner_sha256"], manifest["runner_sha256"])
        self.assertEqual(identity["runner_sha256"], manifest["exporter_sha256"])
        self.assertEqual(str(older_run), manifest["workdir"])

        self.args.genome.write_text(">chr1\nATGAAATAG\n")
        with self.assertRaisesRegex(ValueError, "does not match"):
            run_braker4.export_completed_run(self.args)
        self.assertEqual(manifest, json.loads((self.args.output_dir / "run.json").read_text()))

    def test_bad_results_do_not_replace_previous_exports(self):
        run, identity, cmd, _ = run_braker4.prepare_run(self.args)
        results = self.outputs(run)
        old = self.args.output_dir / "braker.gff3"
        old.write_text("previous annotation\n")
        with gzip.open(results / "braker.aa.gz", "wt") as handle:
            handle.write(">wrong_id\nMK\n")
        with self.assertRaisesRegex(ValueError, "FASTA IDs differ"):
            run_braker4.export_results(run, self.args.output_dir, self.args.sample, identity, cmd)
        self.assertEqual("previous annotation\n", old.read_text())
        self.assertFalse((self.args.output_dir / "run.json").exists())

    def test_missing_collected_file_forces_collection_again(self):
        run, _, _, _ = run_braker4.prepare_run(self.args)
        results = self.outputs(run)
        (results / "braker.aa.gz").unlink()
        resumed, _, command, _ = run_braker4.prepare_run(self.args)
        self.assertEqual(run, resumed)
        self.assertEqual("collect_results", command[command.index("--forcerun") + 1])

    def test_resume_reuses_matching_pinned_predictors(self):
        new_template = self.args.template
        old = configparser.ConfigParser()
        old.read(new_template)
        old.remove_option("PARAMS", "auto_mode_fallback")
        self.args.template = self.root / "old.ini"
        with self.args.template.open("w") as handle:
            old.write(handle)
        self.args.source_manifest.write_text(json.dumps({"commit": run_braker4.RESUME_FROM_COMMIT}))
        run, original, _, _ = run_braker4.prepare_run(self.args)
        sample = self.outputs(run).parent
        for name in ("genome.fa", "augustus.hints.fixed.gtf", "GeneMark-ETP/genemark.gtf", "braker.tsebra.raw.gtf"):
            path = sample / name
            path.parent.mkdir(exist_ok=True)
            path.write_text("saved prediction input\n")
        config_before = (run / "config.ini").read_bytes()
        self.args.template = new_template
        self.args.source_manifest.write_text(json.dumps({"commit": run_braker4.RESUME_TO_COMMIT}))
        self.args.resume_postprocessing_run = run
        actual, identity, command, env = run_braker4.prepare_run(self.args)
        self.assertEqual(actual, run)
        self.assertEqual(original["source"], identity["predictions_source"])
        self.assertEqual(original, json.loads((run / "inputs.json").read_text()))
        self.assertEqual(config_before, (run / "config.ini").read_bytes())
        self.assertNotEqual(env["BRAKER4_CONFIG"], str(run / "config.ini"))
        self.assertEqual("best_by_compleasm", command[command.index("--forcerun") + 1])
        self.assertNotIn("run_genemark_etp", command)
        self.args.genome.write_text(">changed\nACTG\n")
        with self.assertRaisesRegex(ValueError, "matching input files"):
            run_braker4.prepare_run(self.args)

    def test_suppressed_rescue_error_and_missing_comparisons_are_rejected(self):
        run, _, _, _ = run_braker4.prepare_run(self.args)
        result = self.outputs(run)
        braker4_compleasm.validate_rescue(run, self.args.sample)
        log = result.parent / "best_by_compleasm.log"
        log.write_text("STATUS: FAILED: parser exited with code 1\nCompleted all three comparisons.\n")
        with self.assertRaisesRegex(ValueError, "best_by_compleasm did not complete successfully"):
            braker4_compleasm.validate_rescue(run, self.args.sample)
        log.write_text("Completed all three comparisons.\n")
        with self.assertRaisesRegex(ValueError, "best_by_compleasm did not complete successfully"):
            braker4_compleasm.validate_rescue(run, self.args.sample)
        log.write_text("STATUS: OK: original gene set kept\n")
        braker4_compleasm.validate_rescue(run, self.args.sample)
        missing = result.parent / "best_by_compleasm_tmp/single/augustus/summary.txt"
        missing.unlink()
        with self.assertRaisesRegex(ValueError, "Missing Compleasm rescue"):
            braker4_compleasm.validate_rescue(run, self.args.sample)

    def test_compleasm_checks_the_files_extracted_by_the_pinned_downloader(self):
        source = self.root / "lineage"
        candidate = self.root / "private"
        for directory in (source, candidate):
            (directory / "hmms").mkdir(parents=True)
            for name in ("scores_cutoff", "refseq_db.faa.gz", "hmms/test.hmm"):
                (directory / name).write_text(name)
        (source / "dataset.cfg").write_text("name=fixture\n")
        (source / "lengths_cutoff").write_text("test\t100\n")

        braker4_compleasm.verify_lineage(source, candidate)
        (candidate / "scores_cutoff").write_text("changed")
        with self.assertRaisesRegex(ValueError, "scores_cutoff"):
            braker4_compleasm.verify_lineage(source, candidate)

    def test_compleasm_fixed_lineage_catalog_omits_new_placement_names(self):
        md5 = "a" * 32
        lineages = ["eukaryota_odb12", "embryophyta_odb12"]
        rows = [f"{name}\t2026-05-22\t{md5}\tEukaryota\tlineages" for name in lineages]
        placement = f"list_of_reference_markers.eukaryota_odb12.2.txt\t2026-05-22\t{md5}\tEukaryota\tplacement_files"
        catalog = ("\n".join([*rows, placement]) + "\n").encode()
        cache = self.root / "cache"
        with self.assertRaisesRegex(ValueError, "checksum"):
            braker4_compleasm.write_fixed_lineage_catalog(cache, lineages[1], catalog, md5)
        braker4_compleasm.write_fixed_lineage_catalog(cache, lineages[1], catalog, hashlib.md5(catalog).hexdigest())
        self.assertEqual(rows, (cache / "file_versions.tsv").read_text().splitlines())
        self.assertEqual(catalog, (cache / "file_versions.upstream.tsv").read_bytes())

    def test_compleasm_publishes_and_reuses_odb12_cache_without_length_cutoffs(self):
        source = self.args.busco_downloads / "lineages" / self.args.lineage
        (source / "hmms").mkdir()
        entries = ("scores_cutoff", "refseq_db.faa.gz", "hmms/test.hmm")
        for name in entries:
            (source / name).write_text(name)
        catalog = "".join(
            f"{name}\t2026-05-22\t{'a' * 32}\tEukaryota\tlineages\n"
            for name in (self.args.lineage, "eukaryota_odb12")
        ).encode()
        checksum = hashlib.md5(catalog).hexdigest().encode()

        def download(command, **kwargs):
            lineage = Path(command[-1]) / self.args.lineage
            (lineage / "hmms").mkdir(parents=True)
            for name in entries:
                (lineage / name).write_text(name)

        with mock.patch.object(braker4_compleasm, "urlopen", side_effect=[io.BytesIO(checksum), io.BytesIO(catalog)]):
            with mock.patch.object(braker4_compleasm.subprocess, "run", side_effect=download) as runner:
                cache = braker4_compleasm.download_cache(self.args, self.root, ["apptainer"])
                self.assertEqual(cache, braker4_compleasm.download_cache(self.args, self.root, ["apptainer"]))
                self.assertEqual(1, runner.call_count)
        self.assertFalse((cache / self.args.lineage / "lengths_cutoff").exists())

    def test_container_pull_is_atomic_and_reuses_snakemake_image(self):
        config = configparser.ConfigParser()
        config.read(self.args.template)
        uri = config["containers"]["braker3_image"]
        path = self.args.container_cache / (hashlib.md5(uri.encode()).hexdigest() + ".simg")

        def pull(command, **kwargs):
            self.assertEqual(["apptainer", "pull"], command[:2])
            self.assertEqual(uri, command[-1])
            self.assertFalse(path.exists())
            Path(command[2]).write_bytes(b"image")

        with mock.patch.object(braker4_compleasm.subprocess, "run", side_effect=pull) as runner:
            first = braker4_compleasm.container_command(self.args, self.root)
            self.assertEqual(first, braker4_compleasm.container_command(self.args, self.root))
            self.assertEqual(1, runner.call_count)
        self.assertEqual(str(path), first[-2])
        self.assertEqual(b"image", path.read_bytes())

    def test_runner_validates_compleasm_before_export(self):
        run, _, expected, _ = run_braker4.prepare_run(self.args)
        collected = self.outputs(run)
        argv = ["run_braker4"]
        for name in ("source", "source_manifest", "template", "genome", "proteins",
                     "busco_downloads", "container_cache", "output_dir",
                     "sample", "lineage", "threads", "mem_mb"):
            argv.extend(["--" + name.replace("_", "-"), str(getattr(self.args, name))])
        argv.extend(["--rnaseq-r1", *map(str, self.args.rnaseq_r1),
                     "--rnaseq-r2", *map(str, self.args.rnaseq_r2)])
        with mock.patch.object(sys, "argv", argv):
            with mock.patch.object(run_braker4, "prepare_cache") as prepared:
                def execute(command, **kwargs):
                    prepared.assert_called_once()
                    self.assertEqual(expected, command)
                with mock.patch.object(run_braker4.subprocess, "run", side_effect=execute) as child:
                    run_braker4.main()
                    child.assert_called_once()
            manifest = self.args.output_dir / "run.json"
            completed_manifest = manifest.read_bytes()
            (collected.parent / "compleasm_proteins/summary.txt").write_text("Assessment failed\n")
            with mock.patch.object(run_braker4, "prepare_cache"), mock.patch.object(run_braker4.subprocess, "run"):
                with self.assertRaisesRegex(ValueError, "Incomplete Compleasm"):
                    run_braker4.main()
            self.assertEqual(completed_manifest, manifest.read_bytes())

    def test_placeholder_compleasm_summary_is_not_accepted(self):
        path = self.root / "summary.txt"
        path.write_text("Compleasm assessment failed\n")
        with self.assertRaisesRegex(ValueError, "Incomplete Compleasm"):
            braker4_compleasm.complete_summary(path)

    def test_longest_cds_ties_use_input_order(self):
        cds, gff, output = (self.root / name for name in ("cds.fa", "gene.gff3", "longest.fa"))
        cds.write_text(CDS.replace("ATGAAA\n", "ATGAAAAAA\n"))
        gff.write_text(GFF3)
        select_longest_cds.select_longest(cds, gff, output)
        self.assertEqual({"transcriptA", "independent"}, run_braker4.fasta_ids(output))

    def test_archive_checksum_and_extraction(self):
        archive = self.root / "source.tar.gz"
        with tarfile.open(archive, "w:gz") as handle:
            data = b"rule all:\n    input: []\n"
            member = tarfile.TarInfo("BRAKER4-test/Snakefile")
            member.size = len(data)
            handle.addfile(member, io.BytesIO(data))
        destination = self.root / "installed"
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            download_braker4.install_archive(archive, destination, "wrong")
        self.assertFalse(destination.exists())
        digest = download_braker4.sha256_file(archive)
        download_braker4.install_archive(archive, destination, digest)
        self.assertTrue((destination / "Snakefile").exists())

    @unittest.skipUnless(os.environ.get("BRAKER4_TEST_SOURCE"), "Set BRAKER4_TEST_SOURCE to check the real upstream DAG")
    def test_pinned_upstream_etp_dry_run(self):
        self.args.source = Path(os.environ["BRAKER4_TEST_SOURCE"])
        self.args.dry_run = True
        _, _, command, env = run_braker4.prepare_run(self.args)
        result = subprocess.run(command, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(0, result.returncode, result.stdout)
        self.assertIn("run_genemark_etp", result.stdout)
        self.assertIn("convert_gtf_to_gff3", result.stdout)
        self.assertIn("collect_results", result.stdout)

    @unittest.skipUnless(os.environ.get("BRAKER4_TEST_SOURCE"), "Set BRAKER4_TEST_SOURCE to check the real upstream DAG")
    def test_pinned_upstream_varus_etp_dry_run(self):
        self.use_varus()
        self.args.source = Path(os.environ["BRAKER4_TEST_SOURCE"])
        self.args.dry_run = True
        # Exercise argument parsing too: the parent passes empty R1/R2 lists.
        argv = [str(REPO / "workflow/scripts/run_braker4.py")]
        for key, value in vars(self.args).items():
            argv.append("--" + key.replace("_", "-"))
            if isinstance(value, list):
                argv.extend(map(str, value))
            elif value is not True:
                argv.append(str(value))
        result = subprocess.run([sys.executable, *argv], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(0, result.returncode, result.stdout)
        self.assertIn("rule run_varus:", result.stdout)
        self.assertIn("rule run_genemark_etp:", result.stdout)
        self.assertIn("rule collect_results:", result.stdout)
        self.assertNotIn("rule hisat2_align:", result.stdout)
        self.assertNotIn("rule fastp:", result.stdout)


if __name__ == "__main__":
    unittest.main()
