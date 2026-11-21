import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

from check_lmfdb_inessential import (check_field, db_fields, inessential_primes,  # noqa: E402
                                     n_irreducible, parse_primes)
from compare_lmfdb import DiscRange, Log  # noqa: E402
from monogenic import pari  # noqa: E402

C961 = [8, -10, -1, 1]  # 3.3.961.1: 2 splits completely
Q2197 = [3, 4, 2, -1, 1]  # 4.0.2197.1: 3 is inessential
C23 = [-1, -1, 0, 1]  # 3.1.23.1


def kinds(problems):
    return sorted(k for k, _ in problems[0])


class TestInessential(unittest.TestCase):
    def test_n_irreducible(self):
        self.assertEqual([n_irreducible(2, f) for f in range(1, 7)], [2, 1, 2, 3, 6, 9])
        self.assertEqual([n_irreducible(3, f) for f in range(1, 4)], [3, 3, 8])
        for p in (2, 3, 5, 7):
            for f in range(1, 9):
                self.assertEqual(sum(d * n_irreducible(p, d) for d in range(1, f + 1) if f % d == 0), p ** f)

    def test_primes(self):
        self.assertEqual(inessential_primes(C961), [2])
        self.assertEqual(inessential_primes(Q2197), [3])
        self.assertEqual(inessential_primes(C23), [])
        self.assertEqual(inessential_primes([1, -1, 1]), [])  # quadratic: never
        # the compositum of 3.3.961.1 and Q(i): 2 ramifies in Q(i), but there are still 3
        # primes of degree 1 above 2, and only 2 linear polynomials over F_2
        f = pari.polcompositum(pari.Polrev(C961), pari("x^2 + 1"))[0]
        self.assertIn(2, inessential_primes([int(c) for c in pari.Vecrev(f)]))

    def test_agrees_with_monogenic_gp(self):
        # the same criterion is in monogenic.gp (used by --method local); compare on the
        # cubic and quartic fields with |D| <= 2000
        from check_lmfdb_zk import common_index_divisors
        for g in ("[3,1]", "[3,2]", "[4,1]", "[4,2]", "[4,3]", "[4,4]", "[4,5]"):  # C3, S3, C4, ..., S4
            for f in pari.nflist(pari(g), [1, 2000]):
                self.assertEqual(inessential_primes([int(c) for c in pari.Vecrev(f)]),
                                 common_index_divisors(f), str(f))

    def test_parse_primes(self):
        self.assertEqual(parse_primes([3, 2]), [2, 3])
        self.assertEqual(parse_primes([]), [])
        self.assertEqual(parse_primes("[2, 3]"), [2, 3])
        self.assertEqual(parse_primes("{2,3}"), [2, 3])
        self.assertEqual(parse_primes("[]"), [])
        self.assertIsNone(parse_primes(None))
        self.assertIsNone(parse_primes("\\N"))

    def test_check_field(self):
        self.assertEqual(check_field("3.3.961.1", C961, [2], 2, -1, -1), ([], [2]))
        self.assertEqual(check_field("3.1.23.1", C23, [], 1, 1, 1), ([], []))
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [], 2)), ["inessentialp"])
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2, 3], 2)), ["inessentialp"])
        self.assertEqual(kinds(check_field("3.1.23.1", C23, [2])), ["inessentialp"])
        p = check_field("3.3.961.1", C961, None)
        self.assertEqual(kinds(p), ["missing"])
        self.assertIn("it should be {2}", p[0][0][1])
        # the index column: its primes, and it divides [O_K : Z[a]] = 2
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 1)), ["index"])
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 4)), ["index"])
        self.assertEqual(kinds(check_field("4.0.2197.1", Q2197, [3], 3)), [])
        self.assertEqual(kinds(check_field("4.0.2197.1", Q2197, [3], 6)), ["index", "index"])
        # monogenic = 1 is impossible with an inessential prime
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 2, 1, None)), ["monogenic"])
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 2, None, 1)), ["ours"])
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 2, -1, -1)), [])

    def test_db_fields(self):
        fields = [{"label": "3.1.23.1", "coeffs": C23, "inessentialp": [], "index": 1, "monogenic": 1},
                  {"label": "3.3.961.1", "coeffs": C961, "inessentialp": [], "index": 1, "monogenic": 0},
                  {"label": "3.3.473.1", "coeffs": [-1, -10, 1, 1], "monogenic": 0},
                  {"label": "4.0.2197.1", "coeffs": Q2197, "inessentialp": [3], "index": 3, "monogenic": -1}]

        class NF:
            def count(self, query):
                return sum(1 for r in fields if int(r["label"][0]) == query["degree"])

            def search(self, query, cols):
                for r in fields:
                    if int(r["label"][0]) == query["degree"]:
                        yield {k: r[k] for k in cols if k in r}

        db = type("DB", (), {"nf_fields": NF()})()
        rows = list(db_fields(3, DiscRange(), Log(quiet=True), (0, 1), db=db))
        self.assertEqual([r[0] for r in rows], ["3.1.23.1", "3.3.961.1", "3.3.473.1"])
        self.assertEqual([kinds(check_field(*r)) for r in rows], [[], ["index", "inessentialp"], ["missing"]])
        rows = list(db_fields(4, DiscRange(), Log(quiet=True), (0, 1), db=db))
        self.assertEqual([kinds(check_field(*r)) for r in rows], [[]])

    def test_script(self):
        rows = ["label|coeffs|inessentialp|index|monogenic",
                "3.1.23.1|[-1,-1,0,1]|[]|1|1",
                "3.3.961.1|[8,-10,-1,1]|[]|2|\\N",
                "3.3.473.1|[-1,-10,1,1]|\\N",
                "4.0.2197.1|[3,4,2,-1,1]|[3]|3|-1"]
        with tempfile.TemporaryDirectory() as tmp:
            src, out, upd = (os.path.join(tmp, name) for name in ("in.txt", "problems.txt", "update.txt"))
            with open(src, "w") as fh:
                fh.write("\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_lmfdb_inessential.py"),
                                   "--input", src, "--write", out, "--update", upd, "--no-data", "-q"],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn("3.3.961.1: inessentialp is none, but the inessential primes are {2}", proc.stdout)
            self.assertIn("degree 3: 3 fields checked, 1 with inessential primes, 2 with problems; "
                          "1 missing; 1 inessentialp\n", proc.stdout)
            self.assertIn("degree 4: 1 fields checked, 1 with inessential primes, 0 with problems\n", proc.stdout)
            with open(out) as fh:
                self.assertEqual([line.split("|")[:2] for line in fh.read().splitlines()],
                                 [["3.3.473.1", "missing"], ["3.3.961.1", "inessentialp"]])
            with open(upd) as fh:
                self.assertEqual(fh.read(), "label|inessentialp\ntext|integer[]\n\n3.3.473.1|[]\n3.3.961.1|[2]\n")


if __name__ == "__main__":
    unittest.main()
