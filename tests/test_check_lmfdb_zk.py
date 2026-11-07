import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

from check_lmfdb_zk import check_zk, db_fields  # noqa: E402
from compare_lmfdb import DiscRange, Log  # noqa: E402

F = [8673, 0, 252, 0, 1]  # 4.0.19982592.7, x^4 + 252x^2 + 8673
GOOD = ["1", "a", "1/49*a^2 + 18/7", "1/49*a^3 + 25/7*a"]


def kinds(problems):
    return sorted(k for k, _ in problems)


class TestCheckZk(unittest.TestCase):
    def test_good(self):
        self.assertEqual(check_zk("4.0.19982592.7", F, GOOD, 1), [])
        self.assertEqual(check_zk("4.0.19982592.7", F, GOOD, None, nfdisc=True), [])
        self.assertEqual(check_zk("3.1.23.1", [-1, -1, 0, 1], ["1", "a", "a^2"], 1), [])
        # coefficient lists work too, and any basis of O_K passes
        self.assertEqual(check_zk("2.0.3.1", [1, -1, 1], [[1], [0, 1]], 1), [])
        self.assertEqual(check_zk("2.0.3.1", [1, -1, 1], ["1", "a + 5"], 1), [])

    def test_power_basis_is_wrong(self):
        # the LMFDB bug: 1, a, a^2, a^3 is not O_K (index 2401); the index column, the
        # index of the field, is 1 (no common index divisors), and that is right
        p = check_zk("4.0.19982592.7", F, ["1", "a", "a^2", "a^3"], 1)
        self.assertEqual(kinds(p), ["disc"])
        self.assertIn("index 2401 in O_K", dict(p)["disc"])
        # 2401 = 7^4, but 7 is not a common index divisor
        self.assertEqual(kinds(check_zk("4.0.19982592.7", F, GOOD, 2401)), ["index"])

    def test_field_index(self):
        # 2 splits completely in 3.3.961.1 and there are only 2 linear polynomials over F_2,
        # so 2 divides [O_K : Z[alpha]] for every alpha: the index of the field is 2
        f, zk = [8, -10, -1, 1], ["1", "a", "1/2*a^2 - 1/2*a"]
        self.assertEqual(check_zk("3.3.961.1", f, zk, 2), [])
        p = check_zk("3.3.961.1", f, zk, 1)
        self.assertEqual(kinds(p), ["index"])
        self.assertIn("common index divisors (the primes dividing the index of the field) are {2}", p[0][1])
        # 4 has the right primes but does not divide i(a) = 2
        p = check_zk("3.3.961.1", f, zk, 4)
        self.assertEqual(kinds(p), ["index"])
        self.assertIn("does not divide [O_K : Z[a]] = 2", p[0][1])
        # 3 is not a common index divisor, and 6 does not divide 2
        self.assertEqual(kinds(check_zk("3.3.961.1", f, zk, 6)), ["index", "index"])
        # a quartic where 3 is a common index divisor
        self.assertEqual(check_zk("4.0.2197.1", [3, 4, 2, -1, 1], None, 3)[0][0], "format")
        self.assertEqual(kinds(check_zk("4.0.2197.1", [3, 4, 2, -1, 1], None, 1)), ["format", "index"])

    def test_not_integral(self):
        p = check_zk("4.0.1197.2", [9, 3, -5, -1, 1], ["1", "a", "1/2*a^2", "a^3"])
        self.assertEqual(kinds(p), ["disc", "not_integral"])
        self.assertIn("zk element 3 (1/2*a^2) is not", dict(p)["not_integral"])

    def test_format(self):
        self.assertEqual(kinds(check_zk("3.1.23.1", [-1, -1, 0, 1], ["1", "a"])), ["format"])
        self.assertEqual(kinds(check_zk("3.1.23.1", [-1, -1, 0, 1], None)), ["format"])
        self.assertEqual(kinds(check_zk("3.1.23.1", [-1, -1, 0, 1], ["1", "b", "a^2"])), ["format"])
        self.assertEqual(kinds(check_zk("3.1.23.1", [-1, -1, 0, 1], ["1", "a", "2*a"])), ["dependent"])

    def test_wrong_label(self):
        # a basis of O_K for a field whose discriminant is not the label's
        self.assertIn("disc", kinds(check_zk("3.1.31.1", [-1, -1, 0, 1], ["1", "a", "a^2"])))
        self.assertIn("nfdisc", kinds(check_zk("3.1.31.1", [-1, -1, 0, 1], ["1", "a", "a^2"], nfdisc=True)))

    def test_db_fields(self):
        # coeffs and index from nf_fields, zk from nf_fields_extra, joined on the label
        fields = {"4.0.117.1": ([1, 1, -1, -1, 1], 1), "4.0.19982592.7": (F, 1), "3.1.23.1": ([-1, -1, 0, 1], 1)}
        extra = {"4.0.117.1": ["1", "a", "a^2", "a^3"], "4.0.19982592.7": ["1", "a", "a^2", "a^3"]}
        queries = []

        class NF:
            def count(self, query):
                return sum(1 for label in fields if int(label[0]) == query["degree"])

            def search(self, query, cols):
                for label, (coeffs, index) in sorted(fields.items()):
                    if int(label[0]) == query["degree"]:
                        yield {"label": label, "coeffs": coeffs, "index": index}

        class Extra:
            def search(self, query, cols):
                queries.append(query["label"]["$in"])
                return [{"label": label, "zk": extra[label]} for label in query["label"]["$in"] if label in extra]

        db = type("DB", (), {"nf_fields": NF(), "nf_fields_extra": Extra()})()
        rows = list(db_fields(4, DiscRange(), Log(quiet=True), (0, 1), batch=1, db=db))
        self.assertEqual([(r[0], r[2]) for r in rows], [("4.0.117.1", extra["4.0.117.1"]),
                                                       ("4.0.19982592.7", extra["4.0.19982592.7"])])
        self.assertEqual(queries, [["4.0.117.1"], ["4.0.19982592.7"]])
        self.assertEqual([kinds(check_zk(*r)) for r in rows], [[], ["disc"]])

    def test_script(self):
        rows = ['label|coeffs|zk|index',
                '4.0.19982592.7|[8673,0,252,0,1]|["1","a","a^2","a^3"]|1',
                '3.1.23.1|[-1,-1,0,1]|["1","a","a^2"]|\\N']
        with tempfile.TemporaryDirectory() as tmp:
            src, out = os.path.join(tmp, "zk.txt"), os.path.join(tmp, "problems.txt")
            with open(src, "w") as fh:
                fh.write("\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_lmfdb_zk.py"),
                                   "--input", src, "--write", out, "-q"], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn("degree 3: 1 fields checked, 0 with problems", proc.stdout)
            self.assertIn("degree 4: 1 fields checked, 1 with problems; 1 disc\n", proc.stdout)
            with open(out) as fh:
                lines = fh.read().splitlines()
            self.assertEqual([line.split("|")[:2] for line in lines],
                             [["4.0.19982592.7", "disc"]])


if __name__ == "__main__":
    unittest.main()
