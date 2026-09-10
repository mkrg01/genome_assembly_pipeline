"""Regression coverage for the units displayed on linear genome tracks."""

import importlib.util
from pathlib import Path
import runpy
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch


SCRIPTS_DIR = Path(__file__).parents[1] / "workflow" / "scripts"
HAS_PLOTTING = all(
    importlib.util.find_spec(name) is not None for name in ("matplotlib", "pandas")
)


@unittest.skipUnless(HAS_PLOTTING, "Matplotlib and pandas are required")
class LinearPlotTests(unittest.TestCase):
    def test_window_labels_match_each_track_regardless_of_load_order(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contigs = root / "contigs.tsv"
            contigs.write_text("ctg1\t1200000\n")
            gene = root / "gene.bed"
            # A short contig comes first; a terminal partial window comes last.
            gene.write_text(
                "short\t0\t10000\t1\t10\t10000\t0.001\n"
                "ctg1\t0\t500000\t10\t500\t500000\t0.001\n"
                "ctg1\t500000\t1000000\t20\t500\t500000\t0.001\n"
                "ctg1\t1000000\t1200000\t5\t200\t200000\t0.001\n"
            )
            repeat = root / "LTR.bed"
            repeat.write_text(
                "ctg1\t0\t1000000\t30\t1000\t1000000\t0.001\n"
                "ctg1\t1000000\t1200000\t7\t200\t200000\t0.001\n"
            )
            tidk = root / "tidk.tsv"
            tidk.write_text(
                "id\twindow\tforward_repeat_number\treverse_repeat_number\n"
                "ctg1\t10000\t2\t3\n"
                "ctg1\t20000\t1\t2\n"
            )

            class NamedInputs(dict):
                __getattr__ = dict.__getitem__

            inputs = NamedInputs(contig=contigs, gene=gene, LTR=repeat, tidk=tidk)
            expected = {
                "gene": "Count per\n500-kb window",
                "LTR": "Count per\n1-Mb window",
                "tidk": "Count per\n10-kb window",
            }
            for order in (("gene", "LTR", "tidk"), ("tidk", "LTR", "gene")):
                with self.subTest(order=order):
                    job = SimpleNamespace(
                        input=inputs,
                        output=[str(root / "linear.pdf")],
                        config={
                            "circos_plot_tracks": [
                                {"id": name, "label": name, "color": "#4C72B0"}
                                for name in order
                            ]
                        },
                    )
                    script_module = ModuleType("snakemake.script")
                    script_module.snakemake = job
                    figures = []

                    def capture_figure(figure, *args, **kwargs):
                        figures.append(figure)

                    with (
                        plt.rc_context(),
                        patch.dict(sys.modules, {
                            "snakemake": ModuleType("snakemake"),
                            "snakemake.script": script_module,
                        }),
                        patch.object(sys, "path", [str(SCRIPTS_DIR), *sys.path]),
                        patch("matplotlib.figure.Figure.savefig", autospec=True, side_effect=capture_figure),
                    ):
                        runpy.run_path(str(SCRIPTS_DIR / "linear_plot.py"))

                    self.assertEqual(len(figures), 1)
                    labels = [axis.get_ylabel() for axis in figures[0].axes]
                    self.assertEqual(labels, [expected[name] for name in order])


if __name__ == "__main__":
    unittest.main()
