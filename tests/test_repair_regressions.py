import ast
import tempfile
import unittest
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def load_parser():
    tree = ast.parse((ROOT / "bioremediation_gene_miner.py").read_text(encoding="utf-8"))
    nodes = [
        n for n in tree.body
        if (
            isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "DIAMOND_COLS" for t in n.targets)
        ) or (
            isinstance(n, ast.FunctionDef) and n.name == "load_diamond"
        )
    ]
    ns = {"pd": pd, "Path": Path}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "<parser>", "exec"), ns)
    return ns["load_diamond"]

class CoverageTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "hits.tsv"
            path.write_text(text, encoding="utf-8")
            return load_parser()(path)

    def test_gapped_full_reference(self):
        row = self.parse("q\tr\t30\t339\t486\t316\t2e-10\t55\t66.5\t1\t316\n")
        self.assertEqual(row.iloc[0]["reference_coverage_pct"], 100.0)

    def test_partial_reference(self):
        row = self.parse("q\tr\t30\t100\t486\t316\t2e-10\t55\t20\t51\t150\n")
        self.assertEqual(row.iloc[0]["reference_coverage_pct"], 31.6)

    def test_old_format_rejected(self):
        with self.assertRaisesRegex(ValueError, "reference coordinates"):
            self.parse("q\tr\t30\t339\t486\t316\t2e-10\t55\t66.5\n")

    def test_invalid_coordinates_rejected(self):
        with self.assertRaisesRegex(ValueError, "Invalid"):
            self.parse("q\tr\t30\t339\t486\t316\t2e-10\t55\t66.5\t1\t317\n")

    def test_empty_output(self):
        self.assertTrue(self.parse("").empty)

class ReferenceQCTests(unittest.TestCase):
    def check(self, family, protein):
        tree = ast.parse((ROOT / "build_reference_db_v0.2.1.py").read_text(encoding="utf-8"))
        guards = [
            n for n in ast.walk(tree)
            if isinstance(n, ast.If)
            and "QC-REVIEW-LACCASE" in ast.unparse(n)
            and "family_or_target" in ast.unparse(n.test)
        ]
        self.assertEqual(len(guards), 1)
        loop = ast.For(
            target=ast.Name(id="_entry", ctx=ast.Store()),
            iter=ast.List(elts=[ast.Constant(value=1)], ctx=ast.Load()),
            body=[guards[0], ast.parse("accepted.append(acc)").body[0]],
            orelse=[]
        )
        module = ast.fix_missing_locations(ast.Module(body=[loop], type_ignores=[]))
        ns = {
            "acc": "TEST", "e": {}, "query": "test",
            "rule": {"family_or_target": family},
            "protein_name": lambda e: protein,
            "genes": lambda e: "pgeF",
            "qc_review": [], "accepted": [],
            "print": lambda *args: None
        }
        exec(compile(module, "<reference-qc>", "exec"), ns)
        return ns

    def test_pgef_excluded_from_laccase(self):
        result = self.check("laccase", "Peptidoglycan editing factor PgeF")
        self.assertEqual(result["accepted"], [])
        self.assertEqual(len(result["qc_review"]), 1)

    def test_named_laccase_retained(self):
        result = self.check("laccase", "Laccase")
        self.assertEqual(result["accepted"], ["TEST"])

    def test_other_family_unaffected(self):
        result = self.check("flavin reductase", "Flavin reductase")
        self.assertEqual(result["accepted"], ["TEST"])

if __name__ == "__main__":
    unittest.main()
