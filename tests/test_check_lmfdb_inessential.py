import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import cypari2

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

from check_lmfdb_inessential import (check_field, db_fields, field_index,  # noqa: E402
                                     inessential_primes, n_irreducible, parse_primes)
from compare_lmfdb import DiscRange, Log  # noqa: E402
from monogenic import pari  # noqa: E402

C961 = [8, -10, -1, 1]  # 3.3.961.1: 2 splits completely
Q2197 = [3, 4, 2, -1, 1]  # 4.0.2197.1: 3 is inessential
Q11025 = [4, 0, 11, 0, 1]  # 4.0.11025.1: index 4
Q20449 = [36, 0, -1, 0, 1]  # 4.0.20449.1: index 6
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

    def test_field_index(self):
        def index(coeffs, **kw):
            out = 1
            for p, v, lo in field_index(coeffs, **kw):
                self.assertEqual(v, lo)
                out *= p ** v
            return out
        self.assertEqual(index(C961), 2)
        self.assertEqual(index(Q2197), 3)
        self.assertEqual(index(C23), 1)
        self.assertEqual(index(Q11025), 4)
        self.assertEqual(index(Q20449), 6)
        # 2 splits completely in these composita of quadratic fields where 2 splits: 4 and
        # 8 primes of degree 1, so the index is 2^v with v the smallest sum of v(x_i - x_j)
        # over 4 (8) distinct 2-adic integers: 2 (16); and the sextic: 3 + 3 primes, 8
        f4 = pari.polredbest(pari.polcompositum(pari("x^2 - 17"), pari("x^2 - 41"))[0])
        f8 = pari.polredbest(pari.polcompositum(f4, pari("x^2 - 73"))[0])
        f6 = pari.polredbest(pari.polcompositum(pari.Polrev(C961), pari("x^2 - 17"))[0])
        for f, i in ((f4, 2 ** 2), (f8, 2 ** 16 * 3), (f6, 2 ** 8)):
            self.assertEqual(index([int(c) for c in pari.Vecrev(f)]), i)
        # the search over O_K / p^k O_K gives up with bounds lo <= v <= U
        r = pari("mono_field_index_global")(pari.nfinit([f6, [2]]), 2, 1)
        self.assertEqual(r[1], 0)
        self.assertTrue(r[2] <= 8 <= r[0])

    def test_two_searches_agree(self):
        # mono_field_index_global (over O_K / p^k O_K) and mono_field_index_p (split up by
        # the primes above p) are different searches; they agree on the cubic and quartic
        # fields with |D| <= 3000 that have inessential primes, and on some sextics
        g, loc = pari("mono_field_index_global"), pari("mono_field_index_p")
        fields = [f for gr in ("[3,2]", "[4,1]", "[4,2]", "[4,3]", "[4,5]") for f in pari.nflist(pari(gr), [1, 3000])]
        fields += list(pari.nflist(pari("[6,3]"), [1, 10**5]))[:40]
        n = 0
        for f in fields:
            cid, nf = inessential_primes([int(c) for c in pari.Vecrev(f)]), None
            for p in cid:
                nf = nf or pari.nfinit([f, cid])
                a, b = g(nf, p), loc(nf, p)
                self.assertEqual((a[0], a[1], b[1]), (b[0], 1, 1), str(f))
                n += 1
        self.assertGreater(n, 50)

    def test_check_index(self):
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 2)), [])
        p = check_field("3.3.961.1", C961, [2], 4)
        self.assertEqual(kinds(p), ["index"])
        self.assertIn("the index column is 4, but the index of the field is 2", p[0][0][1])
        self.assertEqual(p[2], 2)
        self.assertEqual(kinds(check_field("4.0.11025.1", Q11025, [2], 4)), [])
        self.assertEqual(kinds(check_field("4.0.11025.1", Q11025, [2], 2)), ["index"])
        self.assertEqual(kinds(check_field("4.0.20449.1", Q20449, [2, 3], 6)), [])
        self.assertEqual(kinds(check_field("4.0.20449.1", Q20449, [2, 3], 12)), ["index"])
        # without the exact index, 2 has the right primes and divides [O_K : Z[a]] = 8
        self.assertEqual(kinds(check_field("4.0.11025.1", Q11025, [2], 2, exact_index=False)), [])
        self.assertIsNone(check_field("4.0.11025.1", Q11025, [2], 2, exact_index=False)[2])
        # if the search gives up (here: 1 <= v_2 <= 3), an index within the bounds is undecided
        with mock.patch("check_lmfdb_inessential.field_index", return_value=[(2, 3, 1)]):
            p = check_field("4.0.11025.1", Q11025, [2], 4)
            self.assertEqual(kinds(p), ["index_undecided"])
            self.assertIn("the search gave up with 1 <= v_2 <= 3", p[0][0][1])
            self.assertIsNone(p[2])
            self.assertEqual(kinds(check_field("4.0.11025.1", Q11025, [2], 16)), ["index", "index"])

    def test_parse_primes(self):
        self.assertEqual(parse_primes([3, 2]), [2, 3])
        self.assertEqual(parse_primes([]), [])
        self.assertEqual(parse_primes("[2, 3]"), [2, 3])
        self.assertEqual(parse_primes("{2,3}"), [2, 3])
        self.assertEqual(parse_primes("[]"), [])
        self.assertIsNone(parse_primes(None))
        self.assertIsNone(parse_primes("\\N"))

    def test_check_field(self):
        self.assertEqual(check_field("3.3.961.1", C961, [2], 2, -1, -1), ([], [2], 2))
        self.assertEqual(check_field("3.1.23.1", C23, [], 1, 1, 1), ([], [], 1))
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [], 2)), ["inessentialp"])
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2, 3], 2)), ["inessentialp"])
        self.assertEqual(kinds(check_field("3.1.23.1", C23, [2], 1)), ["inessentialp"])
        p = check_field("3.3.961.1", C961, None, 2)
        self.assertEqual(kinds(p), ["missing"])
        self.assertIn("it should be {2}", p[0][0][1])
        # the index column: its primes, and it divides [O_K : Z[a]] = 2
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 1)), ["index"])
        self.assertEqual(kinds(check_field("3.3.961.1", C961, [2], 4)), ["index"])
        self.assertEqual(kinds(check_field("4.0.2197.1", Q2197, [3], 3)), [])
        self.assertEqual(kinds(check_field("4.0.2197.1", Q2197, [3], 6)), ["index"])
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
        self.assertEqual([kinds(check_field(*r)) for r in rows],
                         [[], ["index", "inessentialp"], ["index_missing", "missing"]])
        rows = list(db_fields(4, DiscRange(), Log(quiet=True), (0, 1), db=db))
        self.assertEqual([kinds(check_field(*r)) for r in rows], [[]])

    def test_script(self):
        rows = ["label|coeffs|inessentialp|index|monogenic",
                "3.1.23.1|[-1,-1,0,1]|[]|1|1",
                "3.3.961.1|[8,-10,-1,1]|[]|2|\\N",
                "3.3.473.1|[-1,-10,1,1]|\\N",
                "4.0.2197.1|[3,4,2,-1,1]|[3]|3|-1"]
        with tempfile.TemporaryDirectory() as tmp:
            src, out, upd, idx = (os.path.join(tmp, name) for name in ("in.txt", "problems.txt", "update.txt", "index.txt"))
            with open(src, "w") as fh:
                fh.write("\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_lmfdb_inessential.py"),
                                   "--input", src, "--write", out, "--update", upd, "--index-update", idx,
                                   "--no-data", "-q"],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn("3.3.961.1: inessentialp is none, but the inessential primes are {2}", proc.stdout)
            self.assertIn("degree 3: 3 fields checked, 1 with inessential primes, 2 with problems; "
                          "1 missing; 1 inessentialp; 1 index_missing\n", proc.stdout)
            self.assertIn("degree 4: 1 fields checked, 1 with inessential primes, 0 with problems\n", proc.stdout)
            with open(out) as fh:
                self.assertEqual([line.split("|")[:2] for line in fh.read().splitlines()],
                                 [["3.3.473.1", "missing"], ["3.3.473.1", "index_missing"],
                                  ["3.3.961.1", "inessentialp"]])
            with open(idx) as fh:
                self.assertEqual(fh.read(), "label|index\ntext|integer\n\n3.3.473.1|1\n")
            with open(upd) as fh:
                self.assertEqual(fh.read(), "label|inessentialp\ntext|integer[]\n\n3.3.473.1|[]\n3.3.961.1|[2]\n")

    def test_cid_from_the_polynomial(self):
        # mono_cid_pol (Dedekind where it applies, nfinit([f, [p]]) elsewhere) agrees with
        # mono_common_index_divisors on nfinit([f, primes < n])
        for gr, B in (("[3,2]", 10**5), ("[4,2]", 10**4), ("[4,3]", 10**4), ("[5,2]", 10**6), ("[6,3]", 10**5)):
            for f in list(pari.nflist(pari(gr), [1, B]))[:150]:
                n = int(pari.poldegree(f))
                nf = pari.nfinit([f, pari.primes([2, n - 1])])
                self.assertEqual(pari("mono_cid_pol")(f), pari("mono_common_index_divisors")(nf), str(f))

    def test_pari_call(self):
        # a stack overflow is retried with a larger maximum stack, up to MAX_STACK
        import monogenic
        size, big = pari.stacksize(), pari.stacksizemax()
        try:
            pari.allocatemem(10**6, 10**6, silent=True)
            with mock.patch.object(monogenic, "MAX_STACK", 2**28):
                self.assertEqual(len(monogenic.pari_call(pari, "vector(2*10^5, i, i^2)")), 2 * 10**5)
            pari.allocatemem(10**6, 10**6, silent=True)
            with mock.patch.object(monogenic, "MAX_STACK", 10**6):
                with self.assertRaises(cypari2.PariError):
                    monogenic.pari_call(pari, "vector(2*10^5, i, i^2)")
        finally:
            pari.allocatemem(size, big, silent=True)

    def test_skip_and_errors(self):
        import check_lmfdb_inessential as cli
        rows = ["label|coeffs|inessentialp|index|monogenic",
                "3.1.23.1|[-1,-1,0,1]|[]|1|1",
                "3.1.31.1|[-1,1,0,1]|[]|1|1",
                "3.3.961.1|[8,-10,-1,1]|[]|2|\\N"]
        with tempfile.TemporaryDirectory() as tmp:
            src, out, upd = (os.path.join(tmp, name) for name in ("in.txt", "problems.txt", "update.txt"))
            with open(src, "w") as fh:
                fh.write("\n".join(rows) + "\n")
            real = cli.check_field

            def boom(label, *args):
                if label == "3.1.31.1":
                    raise cypari2.PariError(28)
                return real(label, *args)
            argv = ["x", "--input", src, "--write", out, "--update", upd, "--no-data", "-q"]
            with mock.patch.object(cli, "check_field", boom), mock.patch.object(sys, "argv", argv), \
                    mock.patch("builtins.print") as printed:
                with self.assertRaises(SystemExit):
                    cli.main()
            text = "\n".join(str(c.args[0]) for c in printed.call_args_list)
            self.assertIn("degree 3: 3 fields checked, 1 with inessential primes, 2 with problems; "
                          "1 inessentialp; 1 error", text)
            with open(out) as fh:
                self.assertEqual([line.split("|")[:2] for line in fh.read().splitlines()],
                                 [["3.1.31.1", "error"], ["3.3.961.1", "inessentialp"]])
            proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_lmfdb_inessential.py"),
                                   "--input", src, "--write", out, "--update", upd, "--no-data", "--skip", "2",
                                   "-q"], capture_output=True, text=True)
            self.assertIn("degree 3: 3 fields checked, 1 with inessential primes, 1 with problems", proc.stdout)
            with open(upd) as fh:  # appended, one header
                self.assertEqual(fh.read(), "label|inessentialp\ntext|integer[]\n\n3.3.961.1|[2]\n3.3.961.1|[2]\n")


if __name__ == "__main__":
    unittest.main()
