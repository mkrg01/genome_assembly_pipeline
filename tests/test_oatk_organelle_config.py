import runpy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


COMMON_RULES = Path(__file__).resolve().parents[1] / "workflow/rules/common.smk"
SCRIPTS = COMMON_RULES.parent.parent / "scripts"


def load_common_rules(**overrides):
    config = {
        "pipeline_version": "v0.0.0",
        "organism_name": "Testus_example",
        "genome_version": "v0.0.0",
        "oatk_lineage": "magnoliopsida",
        **overrides,
    }
    return runpy.run_path(str(COMMON_RULES), init_globals={"config": config})


class OatkOrganelleConfigTest(unittest.TestCase):
    def test_selection_controls_all_organelle_paths(self):
        for selection in (
            ["mitochondrion"],
            ["chloroplast"],
            ["mitochondrion", "chloroplast"],
        ):
            with self.subTest(selection=selection):
                rules = load_common_rules(oatk_organelle=selection)
                self.assertEqual(rules["configured_oatk_organelles"](), selection)
                profiles = rules["oatkdb_path"]()
                outputs = rules["oatk_output_path"]()
                concatemer = rules["concatemer_path"]()
                stats = rules["seqkit_stats_organelle_path"]()
                for organelle, short in (("mitochondrion", "mito"), ("chloroplast", "pltd")):
                    if organelle in selection:
                        self.assertEqual(
                            profiles[f"{short}_fam"],
                            f"results/downloads/oatkdb/magnoliopsida_{short}.fam",
                        )
                        self.assertEqual(
                            outputs[f"{short}_ctg_fasta"],
                            f"results/oatk/oatk/{{assembly_name}}.{short}.ctg.fasta",
                        )
                        self.assertEqual(
                            concatemer[short],
                            f"results/oatk/concatemer/{{assembly_name}}.concatemer.{organelle}.fa",
                        )
                        self.assertEqual(
                            stats[f"{short}_tsv"],
                            f"results/oatk/seqkit/{{assembly_name}}_{organelle}_seqkit_stats.tsv",
                        )
                    else:
                        for paths in (profiles, outputs, concatemer, stats):
                            self.assertFalse(any(short in key for key in paths))
                self.assertIn("utg_final_gfa", outputs)
                self.assertIn("all_organelle", concatemer)

    def test_order_does_not_change_selection_or_output_order(self):
        rules = load_common_rules(oatk_organelle=["chloroplast", "mitochondrion"])
        self.assertEqual(rules["oatk_organelle"], ["mitochondrion", "chloroplast"])
        self.assertEqual(list(rules["oatkdb_path"]()), ["mito_fam", "pltd_fam"])

    def test_omitted_selection_defaults_to_both(self):
        rules = load_common_rules()
        self.assertEqual(rules["configured_oatk_organelles"](), ["mitochondrion", "chloroplast"])

    def test_rejects_invalid_selections(self):
        for value in (
            [],
            ["mitochondrion", "mitochondrion"],
            ["mitochondrion", "mito"],
            ["mito", "pltd"],
            ["nucleus"],
            ["mitochondrion_and_chloroplast"],
            ["mitochondrion", None],
            [["chloroplast"]],
            {"mitochondrion": True},
            True,
            1,
            None,
            "mitochondrion",
            "chloroplast",
            "mitochondrion_and_chloroplast",
            "mito_and_pltd",
            "nucleus",
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "oatk_organelle"):
                    load_common_rules(oatk_organelle=value)


class OatkOrganelleScriptsTest(unittest.TestCase):
    def test_scripts_process_only_selected_organelles(self):
        for selection in (
            ["mitochondrion"],
            ["chloroplast"],
            ["mitochondrion", "chloroplast"],
        ):
            with self.subTest(selection=selection), tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                coverage = root / "coverage.txt"
                coverage.write_text("70\n")
                rules = load_common_rules(oatk_organelle=selection)
                inputs = {**rules["oatkdb_path"](), **rules["oatk_output_path"]()}
                outputs = {**rules["concatemer_path"](), **rules["seqkit_stats_organelle_path"]()}
                job = SimpleNamespace(
                    input=SimpleNamespace(**inputs, minimum_kmer_coverage=coverage, hifi_reads="hifi.fq.gz"),
                    output=SimpleNamespace(
                        **{key: root / key for key in outputs},
                        utg_final_gfa=root / "test.utg.final.gfa",
                    ),
                    params=SimpleNamespace(oatk_organelle=rules["oatk_organelle"]),
                    wildcards=SimpleNamespace(assembly_name="test"),
                    threads=2,
                    log=SimpleNamespace(out=root / "out.log", err=root / "err.log"),
                )
                with mock.patch("subprocess.run", return_value=SimpleNamespace(returncode=0)) as run:
                    runpy.run_path(str(SCRIPTS / "oatk.py"), run_name="__main__", init_globals={"snakemake": job})
                cmd = run.call_args.args[0]
                for organelle, flag, short in (("mitochondrion", "-m", "mito"), ("chloroplast", "-p", "pltd")):
                    self.assertEqual(flag in cmd, organelle in selection)
                    if organelle in selection:
                        self.assertEqual(cmd[cmd.index(flag) + 1], inputs[f"{short}_fam"])

                with mock.patch("subprocess.run") as run:
                    runpy.run_path(str(SCRIPTS / "seqkit_stats_organelle.py"), run_name="__main__", init_globals={"snakemake": job})
                expected_fastas = [inputs[f"{short}_ctg_fasta"] for short in ("mito", "pltd") if f"{short}_fam" in inputs]
                self.assertEqual(
                    [call.args[0][-1] for call in run.call_args_list],
                    [fasta for fasta in expected_fastas for _ in range(2)],
                )

                def write_concat(fasta, prefix, out_fasta, stderr):
                    Path(out_fasta).write_text(f">{prefix}_contig\nACGTACGT\n")

                script = runpy.run_path(str(SCRIPTS / "concatenate_organelle_genome.py"))
                concatenate = script["concatenate_organelle_genome"]
                with mock.patch.dict(concatenate.__globals__, run_concat_replace=mock.Mock(side_effect=write_concat)):
                    concatenate(job)
                self.assertEqual(
                    job.output.all_organelle.read_text(),
                    "".join(f">{organelle}_contig\nACGTACGT\n" for organelle in selection),
                )


if __name__ == "__main__":
    unittest.main()
