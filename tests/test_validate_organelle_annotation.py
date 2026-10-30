import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
import validate_organelle_annotation as script

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False


@unittest.skipUnless(HAVE_BIOPYTHON, "Biopython is required")
class OrganelleValidationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def record(self, sequence="ATGAAATAA", identifier="contig", features=None):
        record = SeqRecord(Seq(sequence), id=identifier, name=identifier, description="")
        record.annotations = {"molecule_type": "DNA", "topology": "circular"}
        if features is None:
            features = [self.cds(SimpleLocation(0, len(sequence), strand=1))]
        record.features = features
        return record

    def cds(self, location, **qualifiers):
        return SeqFeature(location, type="CDS", qualifiers={"gene": ["geneX"], "translation": ["MK"], **qualifiers})

    def validate(self, records=None, fastas=None, decisions=None, **kwargs):
        records = records or [self.record()]
        annotation = self.root/"input.gbk"
        fasta = self.root/"input.fa"
        SeqIO.write(records, annotation, "genbank")
        SeqIO.write(records if fastas is None else fastas, fasta, "fasta")
        path = None
        if decisions is not None:
            path = self.root/"decisions.json"
            path.write_text(json.dumps(decisions))
        before = {p: p.read_bytes() for p in (annotation, fasta, path) if p}
        report = script.run_validation(annotation=annotation, fasta=fasta, decisions=path,
            report_json=self.root/"report.json", report_tsv=self.root/"report.tsv", **kwargs)
        self.assertEqual(before, {p: p.read_bytes() for p in before})
        return report

    def codes(self, report):
        return {i["code"] for i in report["issues"]}

    def test_matching_inputs_and_translation_are_read_only(self):
        report = self.validate()
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["CDS"][0]["translation_status"], "matched")
        script.check_report(self.root/"report.json")

    def test_sequence_record_and_duplicate_ids_are_checked(self):
        for fastas, expected in [([self.record("ATGCCCTAA")], "sequence_mismatch"),
                                  ([self.record(identifier="other")], "record_id_mismatch"),
                                  ([self.record(), self.record()], "duplicate_record_id")]:
            with self.subTest(expected=expected):
                self.assertIn(expected, self.codes(self.validate(fastas=fastas)))

    def test_location_bounds_and_unparsed_locations(self):
        for location, code in [(SimpleLocation(0, 12, strand=1), "location_out_of_bounds"),
                                (SimpleLocation(2, 2, strand=1), "empty_product_location")]:
            with self.subTest(code=code):
                self.assertIn(code, self.codes(self.validate([self.record(features=[self.cds(location)])])))
        self.validate()
        path = self.root/"input.gbk"
        path.write_text(path.read_text().replace("CDS             1..9", "CDS             bad_location"))
        report = script.validate_annotation(annotation=path, fasta=self.root/"input.fa")
        self.assertIn("unparsed_location", self.codes(report))
        self.assertEqual(report["CDS_count"], 1)

    def test_mixed_strand_compound_cds_is_extracted_in_biological_order(self):
        sequence = "ATG" + "C"*7 + "TTT" + "C"*7 + "TAA"
        location = CompoundLocation([SimpleLocation(0,3,strand=1), SimpleLocation(10,13,strand=-1),
                                     SimpleLocation(20,23,strand=1)])
        report = self.validate([self.record(sequence, features=[self.cds(location)])])
        self.assertEqual(report["status"], "pass", report["issues"])
        self.assertEqual(report["CDS"][0]["extracted_length"], 9)

    def test_origin_wrap_and_remote_parts(self):
        origin = CompoundLocation([SimpleLocation(6,9,strand=1),SimpleLocation(0,6,strand=1)])
        self.assertEqual(self.validate([self.record("AAATAAATG", features=[self.cds(origin)])])["status"], "pass")
        remote = CompoundLocation([SimpleLocation(0,3,strand=1),SimpleLocation(0,6,strand=1,ref="other")])
        report = self.validate([self.record("ATG",features=[self.cds(remote)]),self.record("AAATAA",identifier="other",features=[])])
        self.assertEqual(report["status"], "pass", report["issues"])

    def test_genetic_code_and_alternative_initiation(self):
        r = self.record("GTGAAATAA", features=[self.cds(SimpleLocation(0,9,strand=1),transl_table=["11"])])
        self.assertEqual(self.validate([r])["status"], "pass")
        r = self.record("ATGTGATAA", features=[self.cds(SimpleLocation(0,9,strand=1),transl_table=["4"],translation=["MW"])])
        self.assertEqual(self.validate([r])["status"], "pass")
        r.features[0].qualifiers["transl_table"]=["1"]
        self.assertIn("internal_stop", self.codes(self.validate([r])))

    def test_partial_ends_and_codon_start(self):
        r = self.record("AAAATAA", features=[self.cds(SimpleLocation(BeforePosition(0),7,strand=1),codon_start=["2"],translation=["K"])])
        self.assertEqual(self.validate([r])["status"], "pass")
        r = self.record("ATGAAAT", features=[self.cds(SimpleLocation(0,AfterPosition(7),strand=1))])
        self.assertEqual(self.validate([r])["status"], "pass")
        r.features[0].qualifiers["codon_start"]=["4"]
        self.assertIn("invalid_translation_parameters", self.codes(self.validate([r])))

    def test_selenocysteine_and_reverse_strand_transl_except(self):
        for strand in (1,-1):
            seq = Seq("ATGTGATAA")
            if strand == -1:
                seq=seq.reverse_complement()
            location = SimpleLocation(0,9,strand=strand)
            exception="(pos:4..6,aa:Sec)" if strand==1 else "(pos:complement(4..6),aa:Sec)"
            r = self.record(str(seq), features=[self.cds(location,transl_except=[exception],translation=["MU"])])
            self.assertEqual(self.validate([r])["status"], "pass")

    def test_documented_partial_stop_and_invalid_exception(self):
        for ending in ("T","TA"):
            sequence="ATGAAA"+ending
            pos="7" if ending=="T" else "7..8"
            r=self.record(sequence,features=[self.cds(SimpleLocation(0,len(sequence),strand=1),transl_except=[f"(pos:{pos},aa:TERM)"])])
            self.assertEqual(self.validate([r])["status"], "pass")
        r=self.record(features=[self.cds(SimpleLocation(0,9,strand=1),transl_except=["(pos:3..5,aa:Sec)"])])
        self.assertIn("invalid_transl_except",self.codes(self.validate([r])))

    def test_rna_editing_replays_genomic_strand_bases_without_using_feature_index(self):
        for strand in (1,-1):
            seq=Seq("ACGAAATAA")
            if strand==-1:seq=seq.reverse_complement()
            r=self.record(str(seq),features=[self.cds(SimpleLocation(0,9,strand=strand),exception=["RNA editing"])])
            location="1..9" if strand==1 else "complement(1..9)"
            site={"cds_pos":2,"genomic_pos":2 if strand==1 else 8,
                  "genomic_base":"C" if strand==1 else "G","edited_base":"T" if strand==1 else "A"}
            decisions={"feature_decisions":[{"record":"contig","gene":"geneX","location":location,
                       "feature_index":999,"translation_added":True,"accepted_sites":[site]}]}
            report=self.validate([r],decisions=decisions)
            self.assertEqual(report["status"],"pass",report["issues"])
            self.assertEqual(report["CDS"][0]["editing_sites_replayed"],1)
            site["genomic_pos"]=3
            self.assertIn("invalid_editing_site",self.codes(self.validate([r],decisions=decisions)))

    def test_equivalent_compound_location_text_matches_editing_decision(self):
        location=CompoundLocation([SimpleLocation(6,9,strand=-1),SimpleLocation(0,6,strand=-1)])
        r=self.record(str(Seq("ACGAAATAA").reverse_complement()),features=[self.cds(location,exception=["RNA editing"])])
        decisions={"feature_decisions":[{"record":"contig","gene":"geneX",
            "location":"complement(join(1..6,7..9))","translation_added":True,
            "accepted_sites":[{"cds_pos":2,"genomic_pos":8,"genomic_base":"G","edited_base":"A"}]}]}
        self.assertEqual(self.validate([r],decisions=decisions)["status"],"pass")

    def test_pseudogenes_skip_translation_but_still_check_coordinates(self):
        r=self.record("TAGTAGTAG",features=[self.cds(SimpleLocation(0,9,strand=1),pseudo=[""])])
        self.assertEqual(self.validate([r])["status"],"pass")
        r.features[0].location=SimpleLocation(0,12,strand=1)
        self.assertIn("location_out_of_bounds",self.codes(self.validate([r])))

    def test_unreconstructable_exceptions_and_missing_translation_require_review(self):
        for qualifiers in [{"exception":["RNA editing"],"translation":["OTHER"]},
                           {"exception":["reasons given in citation"],"translation":["OTHER"]},
                           {"translation":[]}]:
            with self.subTest(qualifiers=qualifiers):
                r=self.record(features=[self.cds(SimpleLocation(0,9,strand=1),**qualifiers)])
                report=self.validate([r])
                self.assertEqual(report["status"],"review")
                self.assertEqual(report["error_count"],0)
                script.check_report(self.root/"report.json")

    def test_source_only_contig_is_not_a_validation_warning(self):
        r=self.record(features=[SeqFeature(SimpleLocation(0,9),type="source")])
        self.assertEqual(self.validate([r])["status"],"pass")

    def test_error_report_is_preserved_and_release_gate_rejects_it(self):
        r=self.record(features=[self.cds(SimpleLocation(0,9,strand=1),translation=["WRONG"])])
        report=self.validate([r])
        self.assertIn("translation_mismatch",self.codes(report))
        run=subprocess.run([sys.executable,str(SCRIPTS/"validate_organelle_annotation.py"),
                            "--check-report",str(self.root/"report.json")],capture_output=True,text=True)
        self.assertEqual(run.returncode,1)
        self.assertTrue((self.root/"report.tsv").exists())
        self.assertIn("translation_mismatch",(self.root/"report.tsv").read_text())

    def test_stale_validation_report_is_rejected(self):
        self.validate()
        with (self.root/"input.gbk").open('a') as f:f.write('\n')
        with self.assertRaisesRegex(ValueError,"stale"):
            script.check_report(self.root/"report.json")


if __name__ == "__main__":
    unittest.main()
