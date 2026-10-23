"""Check span-only filtering and monochrome self-homology plots at PAF coordinates."""

import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


SCRIPTS_DIR = Path(__file__).parents[1] / "workflow" / "scripts"
HAS_PLOTTING = importlib.util.find_spec("matplotlib") is not None


@unittest.skipUnless(HAS_PLOTTING, "Matplotlib is required")
class SelfDotplotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with patch.object(sys, "path", [str(SCRIPTS_DIR), *sys.path]):
            import plot_self_dotplot
        cls.plot = plot_self_dotplot

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.paf = self.root / "self.paf"
        self.offsets, self.rows, self.total = self.plot.build_contig_offsets(
            [("ctg1", 200_000), ("ctg2", 200_000)]
        )

    @staticmethod
    def record(length=60_000, strand="+", qname="ctg2", qstart=3_000,
               tname="ctg1", tstart=10_000, matches=1_000, block_length=None):
        return "\t".join(map(str, [
            qname, 200_000, qstart, qstart + length, strand,
            tname, 200_000, tstart, tstart + length,
            matches, length if block_length is None else block_length, 0,
        ]))

    def parse(self):
        return self.plot.parse_paf(self.paf, self.offsets, min_alignment_length=10_000)

    def render(self, segments, name):
        from matplotlib.figure import Figure

        figures = []
        savefig = Figure.savefig

        def capture_and_save(figure, *args, **kwargs):
            figure.canvas.draw()
            figures.append(figure)
            return savefig(figure, *args, **kwargs)

        output = self.root / f"{name}.pdf"
        with patch.object(Figure, "savefig", autospec=True, side_effect=capture_and_save):
            self.plot.plot_dotplot(
                self.rows, self.total, segments, output, "Genome (Mb)",
            )
        self.assertTrue(output.read_bytes().startswith(b"%PDF"))
        return figures[0]

    def test_span_cutoff_applies_to_both_sequences_and_not_paf_block_length(self):
        short_query = self.record().split("\t")
        short_query[3] = str(int(short_query[2]) + 49_999)
        short_target = self.record().split("\t")
        short_target[8] = str(int(short_target[7]) + 49_999)
        self.paf.write_text("\n".join([
            self.record(length=49_999, block_length=100_000),
            self.record(length=50_000, block_length=10_000),
            "\t".join(short_query),
            "\t".join(short_target),
        ]) + "\n")
        parsed = self.plot.parse_paf(self.paf, self.offsets, min_alignment_length=50_000)
        self.assertEqual(parsed["total_records"], 4)
        self.assertEqual(parsed["kept_records"], 1)
        self.assertEqual(parsed["filtered_records"]["short_span"], 3)
        self.assertEqual(parsed["segments"], [[(10_000, 203_000), (60_000, 253_000)]])
        self.assertEqual(
            self.plot.parse_paf(self.paf, self.offsets, min_alignment_length=49_999)["kept_records"],
            4,
        )

    def test_identity_tags_and_matching_fraction_do_not_filter_or_recolour(self):
        baseline = None
        for tags in ("", "\tdv:f:0.99", "\tdv:f:nan", "\tNM:i:59000\tcg:Z:60000M"):
            self.paf.write_text(self.record(matches=0) + tags + "\n")
            parsed = self.parse()
            self.assertEqual(parsed["kept_records"], 1)
            if baseline is None:
                baseline = parsed
            self.assertEqual(parsed, baseline)

    def test_forward_and_reverse_coordinates(self):
        self.paf.write_text(self.record() + "\n" + self.record(strand="-") + "\n")
        parsed = self.parse()
        self.assertEqual(parsed["segments"], [
            [(10_000, 203_000), (70_000, 263_000)],
            [(70_000, 203_000), (10_000, 263_000)],
        ])
        figure = self.render(parsed, "monochrome")
        self.assertEqual(len(figure.axes), 1)  # No identity colour bar.
        self.assertFalse(figure.axes[0].get_legend())
        collection = next(item for item in figure.axes[0].collections if item.get_zorder() == 3)
        from matplotlib.colors import to_rgba
        self.assertEqual(collection.get_colors().tolist(), [list(to_rgba("#222222"))])
        self.assertIsNone(collection.get_array())
        self.assertFalse(collection.get_rasterized())
        self.assertNotIn("identity", figure.axes[0].get_title().lower())

    def test_sublinewidth_matches_do_not_appear_perpendicular_to_their_strand(self):
        import numpy as np
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        for strand, direction in (("+", 1), ("-", -1)):
            with self.subTest(strand=strand):
                self.paf.write_text(self.record(strand=strand) + "\n")
                rendered = self.render(self.parse(), f"short_{direction}")
                collection = next(c for c in rendered.axes[0].collections if c.get_zorder() == 3)
                collection.remove()
                # Isolate the production stroke at a known physical length:
                # a 0.2 pt match, shorter than its line width.
                figure = Figure(figsize=(1, 1), dpi=2400)
                FigureCanvasAgg(figure)
                axis = figure.add_axes([0, 0, 1, 1])
                half_span = 0.2 / (72 * 2 * 2**0.5)
                collection.set_segments([[
                    (0.5 - half_span, 0.5 - direction * half_span),
                    (0.5 + half_span, 0.5 + direction * half_span),
                ]])
                collection.set_transform(axis.transData)
                axis.add_collection(collection)
                axis.set(xlim=(0, 1), ylim=(0, 1))
                axis.set_axis_off()
                figure.canvas.draw()
                pixels = np.asarray(figure.canvas.buffer_rgba())[1150:1250, 1150:1250, :3]
                weights = 1 - pixels.mean(axis=2) / 255
                yy, xx = np.indices(weights.shape)
                xx = xx - (xx * weights).sum() / weights.sum()
                yy = -yy + (yy * weights).sum() / weights.sum()
                covariance = (weights * xx * yy).sum() / weights.sum()
                self.assertGreater(direction * covariance, 0)

    def test_both_reciprocal_records_are_drawn_at_their_paf_coordinates(self):
        for strand in ("+", "-"):
            original = self.record(strand=strand).split("\t")
            reciprocal = [*original]
            reciprocal[0:4], reciprocal[5:9] = original[5:9], original[0:4]
            self.paf.write_text("\t".join(original) + "\n" + "\t".join(reciprocal) + "\n")
            parsed = self.parse()
            self.assertEqual(parsed["kept_records"], 2)
            figure = self.render(parsed, f"symmetric_{strand}")
            collection = next(item for item in figure.axes[0].collections if item.get_zorder() == 3)
            first, second = collection.get_segments()
            reflected = first[:, ::-1].tolist()
            self.assertIn(second.tolist(), [reflected, reflected[::-1]])

    def test_self_symmetric_inversion_is_not_drawn_twice(self):
        self.paf.write_text(self.record(strand="-", qname="ctg1", qstart=10_000) + "\n")
        figure = self.render(self.parse(), "inversion")
        collection = next(item for item in figure.axes[0].collections if item.get_zorder() == 3)
        self.assertEqual(len(collection.get_segments()), 1)

    def test_only_exact_self_diagonal_is_removed(self):
        self.paf.write_text("\n".join([
            self.record(qname="ctg1", qstart=10_000),
            self.record(qname="ctg1", qstart=10_700),
        ]) + "\n")
        parsed = self.parse()
        self.assertEqual(parsed["kept_records"], 1)
        self.assertEqual(parsed["filtered_records"]["self_diagonal"], 1)

    def test_empty_inputs_produce_labelled_pdfs(self):
        self.paf.write_text("")
        empty = self.parse()
        figure = self.render(empty, "no_matches")
        self.assertIn("No alignments passed filters", [t.get_text() for t in figure.axes[0].texts])
        self.assertEqual(len(figure.axes), 1)
        self.rows, self.total = [], 0
        figure = self.render(empty, "no_contigs")
        self.assertIn("No long contigs", [t.get_text() for t in figure.axes[0].texts])

    def test_gzipped_paf_is_supported(self):
        self.paf.write_text(self.record() + "\n")
        compressed = self.root / "self.paf.gz"
        with gzip.open(compressed, "wt") as handle:
            handle.write(self.paf.read_text())
        self.assertEqual(
            self.plot.parse_paf(compressed, self.offsets, min_alignment_length=10_000),
            self.parse(),
        )

    def test_invalid_paf_records_are_rejected(self):
        unknown = self.record(qname="absent")
        negative = self.record(qstart=-1)
        for record in ("not\ta\tpaf", unknown, negative, self.record(strand="?")):
            self.paf.write_text(record + "\n")
            with self.subTest(record=record), self.assertRaises(ValueError):
                self.parse()

    def test_cli_writes_monochrome_pdf_and_span_summary_for_untagged_paf(self):
        fasta = self.root / "assembly.fa"
        fasta.write_text(">ctg1\n" + "A" * 200_000 + "\n>ctg2\n" + "C" * 200_000 + "\n")
        self.paf.write_text("\n".join([
            self.record(length=99_999, matches=0),
            self.record(length=100_000, matches=0),
        ]) + "\n")
        pdf, contigs, summary = (self.root / name for name in ("out.pdf", "out.tsv", "out.json"))
        result = subprocess.run([
            sys.executable, str(SCRIPTS_DIR / "plot_self_dotplot.py"),
            "--fasta", str(fasta), "--paf", str(self.paf), "--pdf", str(pdf),
            "--contigs", str(contigs), "--summary", str(summary), "--organism-name", "Testus_example",
        ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(pdf.read_bytes().startswith(b"%PDF"))
        self.assertIn("ctg2\t200000\t200000\t400000", contigs.read_text())
        data = json.loads(summary.read_text())
        self.assertEqual(data["plotted_records"], 1)
        self.assertEqual(data["min_alignment_length"], 100_000)
        self.assertEqual(data["total_records"], 2)
        self.assertEqual(data["filtered_records"]["short_span"], 1)
        self.assertEqual(set(data["filtered_records"]), {"short_span", "self_diagonal"})

    def test_invalid_span_arguments_fail_before_plotting(self):
        args = [
            "--fasta", "in.fa", "--paf", "in.paf", "--pdf", "out.pdf",
            "--contigs", "out.tsv", "--organism-name", "test",
        ]
        for value in ("0", "-1", "nan"):
            with self.subTest(value=value), patch("sys.stderr"):
                with self.assertRaises(SystemExit) as caught:
                    self.plot.parse_args([*args, "--min-alignment-length", value])
                self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
