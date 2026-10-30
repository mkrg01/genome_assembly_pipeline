import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
import organelle_feature_inventory as script

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False


@unittest.skipUnless(HAVE_BIOPYTHON, "Biopython is required")
class FeatureInventoryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def feature(self, gene, start=10, end=100, kind="CDS", **qualifiers):
        return SeqFeature(SimpleLocation(start, end, strand=1), type=kind,
                          qualifiers={"gene": [gene], **qualifiers})

    def write(self, name, *feature_sets):
        records = []
        for i, features in enumerate(feature_sets):
            record = SeqRecord(Seq("A" * 1000), id=f"record{i}", name=f"record{i}")
            record.annotations = {"molecule_type": "DNA", "topology": "circular"}
            record.features = features
            records.append(record)
        path = self.root / name
        SeqIO.write(records, path, "genbank")
        return path

    def rows(self, annotation, *references):
        result = script.build_inventory(annotation, list(references))
        return result, {r["gene"].lower(): r for r in result["inventory"]}

    def test_missing_cds_and_rnas_and_no_gene_cds_double_count(self):
        target = self.write("target.gbk", [self.feature("present", kind="gene"), self.feature("present")])
        reference = self.write("ref.gbk", [self.feature("present"), self.feature("rps12"),
                               self.feature("trnX", kind="tRNA"), self.feature("rrnX", kind="rRNA")])
        result, rows = self.rows(target, reference)
        self.assertEqual(rows["present"]["target_copy_count"], 1)
        self.assertEqual(result["missing_candidate_count"], 3)
        self.assertTrue(all(rows[g]["review_reason"] == "missing_candidate"
                            for g in ("rps12", "trnx", "rrnx")))

    def test_exact_duplicates_are_reported_without_inflating_copy_count(self):
        target = self.write("target.gbk", [self.feature("rbcL", locus_tag=["A"]),
                               self.feature("RBCL", locus_tag=["B"])])
        reference = self.write("ref.gbk", [self.feature("rbcL")])
        result, rows = self.rows(target, reference)
        self.assertEqual(result["duplicate_group_count"], 1)
        self.assertEqual(rows["rbcl"]["target_copy_count"], 1)
        self.assertEqual(rows["rbcl"]["target_feature_count"], 2)
        self.assertEqual(rows["rbcl"]["review_reason"], "duplicate_annotation")
        self.assertEqual(result["duplicates"][0]["locus_tags"], ["A", "B"])

    def test_distinct_ir_copies_and_shared_trans_splicing_exon_are_not_duplicates(self):
        def rps12(other_start):
            return SeqFeature(CompoundLocation([
                SimpleLocation(10, 40, strand=1),
                SimpleLocation(other_start, other_start+30, strand=-1),
            ]), type="CDS", qualifiers={"gene": ["rps12"], "trans_splicing": [""]})
        fs = [rps12(200), rps12(500), self.feature("rps12", 10, 40, kind="exon")]
        target = self.write("target.gbk", fs)
        reference = self.write("ref.gbk", fs)
        result, rows = self.rows(target, reference)
        self.assertEqual(rows["rps12"]["target_copy_count"], 2)
        self.assertEqual(result["duplicate_group_count"], 0)
        self.assertFalse(rows["rps12"]["manual_review_candidate"])

    def test_identical_coordinates_on_separate_contigs_are_distinct_copies(self):
        target = self.write("target.gbk", [self.feature("geneX")], [self.feature("geneX")])
        reference = self.write("ref.gbk", [self.feature("geneX")])
        result, rows = self.rows(target, reference)
        self.assertEqual(result["duplicate_group_count"], 0)
        self.assertEqual(rows["genex"]["target_copy_count"], 2)
        self.assertEqual(rows["genex"]["review_reason"], "copy_number_difference")

    def test_reference_files_are_independent_genomes_with_an_observed_count_range(self):
        ref1 = self.write("ref1.gbk", [self.feature("geneX"), self.feature("variable")])
        ref2 = self.write("ref2.gbk", [self.feature("geneX", 10+i*200, 100+i*200) for i in range(3)])
        target = self.write("target.gbk", [self.feature("geneX", 10+i*200, 100+i*200) for i in range(2)])
        _, rows = self.rows(target, ref1, ref2)
        self.assertFalse(rows["genex"]["manual_review_candidate"])
        self.assertEqual(rows["variable"]["reference_present_count"], 1)
        self.assertEqual(rows["variable"]["review_reason"], "missing_candidate")

    def test_no_references_still_checks_duplicates_without_claiming_gene_absence(self):
        target = self.write("target.gbk", [self.feature("geneX"), self.feature("geneX")])
        result, rows = self.rows(target)
        self.assertFalse(result["reference_comparison_evaluated"])
        self.assertIsNone(rows["genex"]["reference_copy_min"])
        self.assertEqual(rows["genex"]["reference_comparison"], "not_evaluated")
        self.assertEqual(rows["genex"]["review_reason"], "duplicate_annotation")

    def test_gene_only_pseudogene_is_present_and_parent_gene_is_not_counted_twice(self):
        target = self.write("target.gbk", [self.feature("geneX", kind="gene", pseudo=[""]),
                               self.feature("other", kind="gene", pseudo=[""]),
                               self.feature("other", pseudo=[""])])
        reference = self.write("ref.gbk", [self.feature("geneX"), self.feature("other")])
        result, rows = self.rows(target, reference)
        self.assertEqual(result["missing_candidate_count"], 0)
        self.assertEqual(rows["genex"]["target_pseudogene_count"], 1)
        self.assertEqual(rows["other"]["target_copy_count"], 1)

    def test_target_only_and_unnamed_features_are_not_called_copy_number_errors(self):
        unknown = self.feature("unused", kind="tRNA")
        del unknown.qualifiers["gene"]
        target = self.write("target.gbk", [self.feature("novel"), unknown])
        reference = self.write("ref.gbk", [self.feature("reference")])
        result, rows = self.rows(target, reference)
        self.assertEqual(rows["novel"]["reference_comparison"], "target_only")
        self.assertFalse(rows["novel"]["manual_review_candidate"])
        self.assertEqual(result["unnamed_target_product_features"], 1)

    def test_reports_are_read_only_and_no_reference_is_explicit(self):
        target = self.write("target.gbk", [self.feature("geneX")])
        original = target.read_bytes()
        result = script.run_inventory(annotation=target, reference_dir=None,
            inventory_tsv=self.root/"inventory.tsv", duplicates_tsv=self.root/"duplicates.tsv",
            summary_json=self.root/"summary.json")
        self.assertEqual(target.read_bytes(), original)
        self.assertEqual(json.loads((self.root/"summary.json").read_text()), result)
        self.assertIn("target_copy_count", (self.root/"inventory.tsv").read_text())
        (self.root/"empty").mkdir()
        with self.assertRaisesRegex(ValueError, "No GenBank reference"):
            script.reference_files(self.root/"empty")


if __name__ == "__main__":
    unittest.main()
