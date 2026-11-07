import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import monogenic  # noqa: E402
from monogenic import (DATA_DIR, cubic_generator, factor_generator, find_generator, format_coeffs,  # noqa: E402
                       label_disc,
                       merge_data, obstruction, parse_coeffs, pari, read_data, root_generates,
                       verify_generator, write_data)

# Examples from LMFDB issue 6783: monogenic, but the root of the polredabs
# polynomial does not generate O_K.
ISSUE_EXAMPLES = {
    "3.1.244.1": [-6, 1, 0, 1],
    "3.1.356.1": [7, 1, -1, 1],
    "3.1.424.1": [8, 0, -1, 1],
    "4.0.1040.2": [4, -4, 1, 0, 1],
    "4.0.1413.1": [9, -3, 7, -1, 1],
    "4.2.2092.1": [-2, -3, 1, -1, 1],
}


class TestGenerators(unittest.TestCase):
    def test_issue_examples(self):
        for label, coeffs in ISSUE_EXAMPLES.items():
            with self.subTest(label=label):
                f = pari.Polrev(coeffs)
                self.assertEqual(pari.nfdisc(f), label_disc(label))
                self.assertNotEqual(pari.poldisc(f), label_disc(label))  # root is not a generator
                g = find_generator(coeffs)
                self.assertIsNotNone(g)
                self.assertEqual(pari.poldisc(pari.Polrev(g)), label_disc(label))
                self.assertTrue(pari.nfisisom(f, pari.Polrev(g)))

    def test_known_generators(self):
        # Pinned so that a change to the canonical choice is noticed.
        self.assertEqual(find_generator(ISSUE_EXAMPLES["3.1.244.1"]), [6, -4, -1, 1])
        # [10,-10,3,-2,1] has the same T2; S(P) = (2,-2,3,3,6,6,2,2) comes first
        self.assertEqual(find_generator(ISSUE_EXAMPLES["4.0.1040.2"]), [2, 6, 3, -2, 1])

    def test_key(self):
        key = lambda g: [int(c) for c in monogenic.generator_key(g)]
        # T2 counts both complex roots: T2(x^2 + 1) = 2, and S(P) = (|a1|, a1, |a2|, a2)
        self.assertEqual(key([1, 0, 1]), [2000000, 0, 0, 1, 1])
        self.assertEqual(key([-2, 3, -1, 1]), key([-2, 3, -1, 1])[:1] + [1, -1, 3, 3, 2, -2])
        # the polredabs order: a1 = -1 comes before a1 = 1 (S = (1,-1,...) < (1,1,...))
        self.assertTrue(monogenic.is_smaller(monogenic.generator_key([-1, -1, 1]),
                                             monogenic.generator_key([-1, 1, 1])))

    def test_search_keeps_root(self):
        # Z[2^(1/3)] is the full ring of integers, so the polynomial itself is kept.
        self.assertEqual(find_generator([-2, 0, 0, 1]), [-2, 0, 0, 1])

    def test_root_generates(self):
        self.assertTrue(root_generates([-2, 0, 0, 1], "3.1.108.1"))
        self.assertTrue(root_generates([1, -1, 1], "2.0.3.1"))
        for label, coeffs in ISSUE_EXAMPLES.items():
            with self.subTest(label=label):
                self.assertFalse(root_generates(coeffs, label))

    def test_dedekind_non_monogenic(self):
        # Dedekind's cubic x^3 - x^2 - 2x - 8 is not monogenic (2 splits completely).
        self.assertIsNone(find_generator([-8, -2, -1, 1]))


class TestNonMonogenic(unittest.TestCase):
    def test_common_index_divisor(self):
        self.assertEqual(obstruction([-8, -2, -1, 1]), "cid 2")  # Dedekind's cubic
        self.assertEqual(obstruction([2, 1, -1, 2, -1, 1]), "cid 2")  # quintic, disc 17161

    def test_local(self):
        # Q(20^(1/3)): the index form 5x^3 - 2y^3 is never +-1 mod 9
        self.assertEqual(obstruction([-20, 0, 0, 1]), "local 9 9")
        # quintic field of discriminant 140625 = 3^2 5^6: nothing is +-1 mod 25
        self.assertEqual(obstruction([-3, 0, -5, 0, 0, 1]), "local 25 25")

    def test_monogenic_fields_have_no_proof(self):
        for label, coeffs in ISSUE_EXAMPLES.items():
            with self.subTest(label=label):
                self.assertIsNone(obstruction(coeffs))
        self.assertIsNone(obstruction([-2, 0, 0, 1]))

    def test_cubic_thue(self):
        for label, coeffs in ISSUE_EXAMPLES.items():
            if label.startswith("3."):
                with self.subTest(label=label):
                    self.assertEqual(cubic_generator(coeffs), find_generator(coeffs))
        self.assertEqual(cubic_generator([-2, 0, 0, 1]), [-2, 0, 0, 1])
        self.assertIsNone(cubic_generator([-8, -2, -1, 1]))
        self.assertIsNone(cubic_generator([-20, 0, 0, 1]))

    def test_cubic_fields_agree(self):
        # Every cubic field with |D| <= 2000: a proof of non-monogenicity must agree
        # with the exact Thue answer, and the search must not beat it.
        for f in pari('concat(nflist("S3", [1, 2000]), nflist("C3", [1, 2000]))'):
            coeffs = [int(c) for c in pari.Vecrev(pari.polredabs(f))]
            gen = cubic_generator(coeffs)
            with self.subTest(coeffs=coeffs):
                if obstruction(coeffs) is not None:
                    self.assertIsNone(gen)
                if gen is None:
                    self.assertIsNone(find_generator(coeffs, 2))


class TestFactorIndexForm(unittest.TestCase):
    def test_examples(self):
        # 4.0.1197.2: no common index divisor or local obstruction, but the index form is
        # (x3^2 - x3 x4 + x4^2) * (a quartic), and that quartic is never +-1 on the six
        # points where the first factor is 1.
        self.assertEqual(factor_generator([9, 3, -5, -1, 1]), (True, None))
        self.assertIsNone(obstruction([9, 3, -5, -1, 1]))
        self.assertEqual(factor_generator([7, 5, -4, -2, 1]), (True, [7, 5, -4, -2, 1]))  # 4.0.1197.1, the root
        # issue examples, with the generators Magma gives (stored in data/degree4.txt)
        self.assertEqual(factor_generator(ISSUE_EXAMPLES["4.0.1040.2"]), (True, [2, 6, 3, -2, 1]))
        self.assertEqual(factor_generator(ISSUE_EXAMPLES["4.0.1413.1"]), (True, [13, -9, -2, 0, 1]))

    def test_real_subfield(self):
        # D4 fields with a real quadratic subfield: no definite factor, but the line and
        # conic in the (U, V) plane (cubic resolvent) have no admissible common point
        self.assertEqual(factor_generator([-4, 0, -2, 0, 1]), (True, None))  # 4.2.1600.1
        self.assertIsNone(obstruction([-4, 0, -2, 0, 1]))
        self.assertEqual(factor_generator([4, 2, -5, -1, 1]), (True, None))  # 4.4.2225.1
        self.assertEqual(factor_generator([-1, -3, 0, -1, 1]), (True, None))  # 4.2.775.1
        # here a candidate (U, V) is left, and U(x), V(x) never take it mod 16
        self.assertEqual(factor_generator([6, -6, -9, 0, 1]), (True, None))

    def test_totally_complex(self):
        # S4: the cubic resolvent is irreducible, so a Thue equation gives (U, V); the
        # field is totally complex, so a definite U - sV leaves finitely many x
        self.assertEqual(factor_generator([4, 1, -1, -1, 1]), (True, None))  # 4.0.3636.1
        self.assertIsNone(obstruction([4, 1, -1, -1, 1]))
        self.assertIsNone(find_generator([4, 1, -1, -1, 1], 3))
        # generators stored from Magma
        self.assertEqual(factor_generator([4, -4, 2, -1, 1]), (True, [8, 3, -5, -1, 1]))  # 4.0.2292.1
        self.assertEqual(factor_generator([4, -2, -1, 0, 1]), (True, [7, -3, 6, -1, 1]))  # 4.0.2920.1

    def test_conic(self):
        # Not totally complex: each (U, V) left by the Thue equation is solved exactly
        # through the conic v U - u V = 0 and a quartic Thue equation.
        self.assertEqual(factor_generator([-9, 1, 3, -1, 1]), (True, None))  # 4.2.10091.1, S4
        self.assertIsNone(obstruction([-9, 1, 3, -1, 1]))
        # generators stored from Magma
        self.assertEqual(factor_generator([-12, 16, -5, -1, 1]), (True, [18, 10, -5, -1, 1]))  # 4.2.28856.1
        self.assertEqual(factor_generator([28, 12, -11, -2, 1]), (True, [-2062, -1132, -179, -2, 1]))  # 4.4.31808.1

    def test_not_applicable(self):
        self.assertEqual(factor_generator([-2, 0, 0, 1]), (True, [-2, 0, 0, 1]))  # root generates
        self.assertEqual(factor_generator([-20, 0, 0, 1]), (False, None))  # cubic

    def test_quartic_fields_agree(self):
        # Totally complex quartic fields with a quadratic subfield, |D| <= 1500: the
        # answer must agree with the proofs of non-monogenicity and the search.
        for f in pari('concat([nflist("C4", [1, 1500], 2), nflist("V4", [1, 1500], 2), nflist("D4", [1, 1500], 2)])'):
            coeffs = [int(c) for c in pari.Vecrev(pari.polredabs(f))]
            decided, gen = factor_generator(coeffs)
            with self.subTest(coeffs=coeffs):
                if not decided:
                    continue
                if gen is None:
                    self.assertIsNone(find_generator(coeffs, 2))
                else:
                    self.assertIsNone(obstruction(coeffs))
                    self.assertTrue(verify_generator(gen, "4.0.%d.1" % abs(int(pari.nfdisc(f)))))
                    found = find_generator(coeffs, 2)
                    if found is not None:
                        self.assertFalse(monogenic.more_canonical(found, gen))


class TestCheckData(unittest.TestCase):
    def test_check_generator(self):
        from check_data import check_generator
        self.assertEqual(check_generator("3.1.244.1", 3, [6, -4, -1, 1], [-6, 1, 0, 1]), [])
        self.assertIn("index 2", check_generator("3.1.244.1", 3, [-6, 1, 0, 1], None)[0])
        # right discriminant, wrong field
        self.assertIn("different field",
                      check_generator("4.0.1040.1", 4, [10, -10, 3, -2, 1], [5, -2, -3, 0, 1])[0])
        # the root of the LMFDB polynomial generates O_K, so it must be stored
        self.assertIn("should be [1,0,-1,1]", check_generator("3.1.23.1", 3, [-1, -1, 0, 1], [1, 0, -1, 1])[0])

    def test_cross_check(self):
        from check_data import cross_check
        self.assertEqual(cross_check("3.1.244.1", 1, [6, -4, -1, 1], [-6, 1, 0, 1], 2), [])
        self.assertIn("marked not monogenic", cross_check("3.1.244.1", -1, None, [-6, 1, 0, 1], 2)[0])
        self.assertEqual(cross_check("3.1.4.1", -1, None, [-8, -2, -1, 1], 3), [])  # Dedekind

    def test_reprove(self):
        from check_data import reprove
        self.assertEqual(reprove("3.1.244.1", 3, 1, [6, -4, -1, 1], [-6, 1, 0, 1]), ([], True))
        self.assertIn("marked not monogenic", reprove("3.1.244.1", 3, -1, None, [-6, 1, 0, 1])[0][0])
        self.assertIn("not the stored one", reprove("3.1.244.1", 3, 1, [-6, 1, 0, 1], [-6, 1, 0, 1])[0][0])
        self.assertIn("no solutions", reprove("3.3.1957.1", 3, 1, [-8, -2, -1, 1], [-8, -2, -1, 1])[0][0])
        self.assertEqual(reprove("5.1.140625.1", 5, -1, None, [-3, 0, -5, 0, 0, 1]), ([], True))
        self.assertEqual(reprove("4.0.1197.2", 4, -1, None, [9, 3, -5, -1, 1]), ([], True))
        self.assertIn("no solutions", reprove("4.0.1197.2", 4, 1, [9, 3, -5, -1, 1], [9, 3, -5, -1, 1])[0][0])
        self.assertEqual(reprove("4.0.1040.2", 4, 1, [2, 6, 3, -2, 1], [4, -4, 1, 0, 1]), ([], True))


class TestMerge(unittest.TestCase):
    """A verified generator always beats -1, and every 1 vs -1 conflict is logged."""
    GEN = [6, -4, -1, 1]   # generates O_K for 3.1.244.1
    BAD = [-6, 1, 0, 1]    # the polredabs polynomial of 3.1.244.1: index 2

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "degree3.txt")
        self.log = os.path.join(self.tmp.name, "conflicts.txt")

    def tearDown(self):
        self.tmp.cleanup()

    def merge(self, stored, new, replace=False):
        write_data(self.path, [("3.1.244.1",) + stored])
        with contextlib.redirect_stderr(io.StringIO()) as err:
            counts = merge_data(self.path, {"3.1.244.1": ("3.1.244.1",) + new}, replace,
                                source="test", log=self.log)
        return read_data(self.path)[0][1:], counts, err.getvalue()

    def logged(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return fh.read().splitlines()

    def test_verify_generator(self):
        self.assertTrue(verify_generator(self.GEN, "3.1.244.1"))
        self.assertFalse(verify_generator(self.BAD, "3.1.244.1"))
        self.assertFalse(verify_generator([-2, 0, 0, 1], "3.1.244.1"))  # another field
        self.assertFalse(verify_generator([-2, 0, 0, 2], "3.1.108.1"))  # not monic
        self.assertFalse(verify_generator(None, "3.1.244.1"))

    def test_generator_overwrites_false_minus_one(self):
        for replace in (False, True):
            with self.subTest(replace=replace):
                row, counts, err = self.merge((-1, None), (1, self.GEN), replace)
                self.assertEqual(row, (1, self.GEN))
                self.assertEqual(counts, (0, 1, 0, 1))
                self.assertIn("CONFLICT 3.1.244.1", err)
        self.assertEqual(len(self.logged()), 2)
        self.assertIn("has -1, the new row (from test) has 1 [6,-4,-1,1]", self.logged()[0])
        self.assertTrue(self.logged()[0].endswith("so stored 1"))

    def test_minus_one_never_overwrites_generator(self):
        # e.g. Magma, which runs with replace=True
        row, counts, _ = self.merge((1, self.GEN), (-1, None), replace=True)
        self.assertEqual(row, (1, self.GEN))
        self.assertEqual(counts, (0, 0, 1, 1))
        self.assertEqual(len(self.logged()), 1)

    def test_bad_generator_loses(self):
        row, counts, _ = self.merge((-1, None), (1, self.BAD))
        self.assertEqual(row, (-1, None))
        self.assertEqual(counts, (0, 0, 1, 1))
        self.assertIn("FAILED verification, so stored -1", self.logged()[0])
        row, counts, _ = self.merge((1, self.BAD), (-1, None))
        self.assertEqual(row, (-1, None))
        self.assertEqual(counts, (0, 1, 0, 1))

    def test_no_conflict(self):
        self.assertEqual(self.merge((0, None), (1, self.GEN))[:2], ((1, self.GEN), (0, 1, 0, 0)))
        self.assertEqual(self.merge((0, None), (-1, None))[:2], ((-1, None), (0, 1, 0, 0)))
        self.assertEqual(self.merge((1, self.GEN), (0, None), True)[:2], ((1, self.GEN), (0, 0, 1, 0)))
        self.assertEqual(self.merge((-1, None), (-1, None))[:2], ((-1, None), (0, 0, 1, 0)))
        self.assertEqual(self.logged(), [])

    def test_write_data_is_atomic(self):
        # A write that fails part way (here: a bad label, found while sorting) leaves the
        # old file as it was, and no temporary file behind.
        write_data(self.path, [("3.1.23.1", 1, [1, 0, -1, 1])])
        with open(self.path) as fh:
            before = fh.read()
        with self.assertRaises(ValueError):
            write_data(self.path, [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.oops.1", -1, None)])
        with open(self.path) as fh:
            self.assertEqual(fh.read(), before)
        self.assertEqual(sorted(os.listdir(self.tmp.name)), ["degree3.txt"])

    def test_parallel_merges_keep_every_row(self):
        # Four processes merge 25 batches of rows each into the same file at the same
        # time, as find_generators.py --shard i/4 does. With the lock, no row is lost.
        script = """
import sys
sys.path.insert(0, %r)
from monogenic import merge_data
path, k = sys.argv[1], int(sys.argv[2])
for b in range(25):
    label = "3.1.%%d.1" %% (1000 * k + b + 1)
    merge_data(path, {label: (label, -1, None)}, source="test", log=%r)
""" % (os.path.dirname(monogenic.__file__), self.log)
        write_data(self.path, [])
        procs = [subprocess.Popen([sys.executable, "-c", script, self.path, str(k)]) for k in range(4)]
        self.assertEqual([p.wait() for p in procs], [0] * 4)
        rows = read_data(self.path)
        self.assertEqual(len(rows), 100)
        self.assertEqual(sorted(row[0] for row in rows),
                         sorted("3.1.%d.1" % (1000 * k + b + 1) for k in range(4) for b in range(25)))

    def test_generate_quadratic_keeps_rows(self):
        # generate_quadratic.py merges: rows outside its range stay, rows in it are added.
        import generate_quadratic
        path = os.path.join(self.tmp.name, "degree2.txt")
        big = ("2.2.1000005.1", 1, [-250001, -1, 1])
        write_data(path, [("2.0.3.1", 1, [1, -1, 1]), big])
        saved = sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG
        try:
            sys.argv = ["generate_quadratic.py", "--max-disc", "20"]
            monogenic.DATA_DIR, monogenic.CONFLICT_LOG = self.tmp.name, self.log
            with contextlib.redirect_stdout(io.StringIO()) as out:
                generate_quadratic.main()
        finally:
            sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = saved
        rows = read_data(path)
        self.assertIn(big, rows)
        self.assertEqual([row[0] for row in rows],
                         ["2.0.3.1", "2.0.4.1", "2.2.5.1", "2.0.7.1", "2.0.8.1", "2.2.8.1", "2.0.11.1",
                          "2.2.12.1", "2.2.13.1", "2.0.15.1", "2.2.17.1", "2.0.19.1", "2.0.20.1",
                          "2.2.1000005.1"])
        self.assertIn("13 quadratic fields: 12 added, 0 replaced, 1 kept", out.getvalue())
        self.assertEqual(self.logged(), [])

    def merge23(self, stored, new, replace=False):
        """Merge generator `new` for 3.1.23.1 into a file storing `stored`."""
        write_data(self.path, [("3.1.23.1", 1, stored)])
        with contextlib.redirect_stderr(io.StringIO()) as err:
            counts = merge_data(self.path, {"3.1.23.1": ("3.1.23.1", 1, new)}, replace,
                                source="test", log=self.log)
        return read_data(self.path)[0][2], counts, err.getvalue()

    def test_canonical_generator_wins(self):
        # Generators of O_K for 3.1.23.1, in the canonical order: the polredabs root,
        # then by T2 norm (3.22, 3.26, 3.83).
        root, g1, g2 = [1, 0, -1, 1], [-1, -1, 0, 1], [-1, 2, -1, 1]
        self.assertTrue(monogenic.more_canonical(root, g1))
        self.assertTrue(monogenic.more_canonical(g1, g2))
        self.assertFalse(monogenic.more_canonical(g2, g1))
        # A smaller generator replaces a larger one, whichever method finds it ...
        coeffs, counts, err = self.merge23(g2, g1)
        self.assertEqual((coeffs, counts), (g1, (0, 1, 0, 0)))
        self.assertIn("BETTER 3.1.23.1: replaced [-1,2,-1,1] by [-1,-1,0,1]", err)
        # ... and the polredabs root replaces any other ...
        self.assertEqual(self.merge23(g1, root)[:2], (root, (0, 1, 0, 0)))
        # ... but nothing replaces a smaller one, with or without replace = True
        # (unless the new generator also checks out).
        self.assertEqual(self.merge23(g1, g2)[:2], (g1, (0, 0, 1, 0)))
        self.assertEqual(self.merge23(root, g1)[:2], (root, (0, 0, 1, 0)))
        self.assertEqual(self.merge23(root, g1, True)[:2], (g1, (0, 1, 0, 0)))
        self.assertEqual(self.logged(), [])

    def test_bad_smaller_generator_rejected(self):
        # The polredabs root of 3.1.244.1 comes first in the canonical order, but Z[alpha]
        # has index 2, so it must not replace the stored generator, even with replace.
        write_data(self.path, [("3.1.244.1", 1, self.GEN)])
        with contextlib.redirect_stderr(io.StringIO()) as err:
            counts = merge_data(self.path, {"3.1.244.1": ("3.1.244.1", 1, self.BAD)}, True,
                                source="test", log=self.log)
        self.assertEqual(read_data(self.path)[0][2], self.GEN)
        self.assertEqual(counts, (0, 0, 1, 1))
        self.assertIn("FAILED verification, so kept the stored one", self.logged()[0])

    def test_recheck_negative(self):
        # A false -1 (as from a buggy Magma run) is corrected by the PARI search, and a
        # true -1 (Dedekind's cubic field) is left alone, as is the 1 row.
        import find_generators
        write_data(self.path, [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.244.1", -1, None),
                               ("3.1.503.1", -1, None)])
        fields = os.path.join(self.tmp.name, "fields.txt")
        with open(fields, "w") as fh:
            fh.write("3.1.23.1|[1,0,-1,1]\n3.1.244.1|[-6,1,0,1]\n3.1.503.1|[8,2,-1,1]\n")
        argv = ["find_generators.py", "--degrees", "3", "--method", "pari", "--recheck-negative",
                "--input", fields]
        saved = sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG
        try:
            sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = argv, self.tmp.name, self.log
            with contextlib.redirect_stdout(io.StringIO()) as out, \
                    contextlib.redirect_stderr(io.StringIO()):
                find_generators.main()
        finally:
            sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = saved
        self.assertEqual(read_data(self.path), [("3.1.23.1", 1, [1, 0, -1, 1]),
                                                ("3.1.244.1", 1, self.GEN), ("3.1.503.1", -1, None)])
        self.assertIn("2 fields checked", out.getvalue())
        self.assertIn("1 CONFLICTS", out.getvalue())
        self.assertEqual(len(self.logged()), 1)
        self.assertIn("find_generators.py --degrees 3 --method pari --recheck-negative", self.logged()[0])


    def test_recheck_generators(self):
        # A stored generator that is not the canonical one is replaced; the canonical
        # one for 3.1.244.1 is left alone.
        import find_generators
        write_data(self.path, [("3.1.23.1", 1, [-1, 2, -1, 1]), ("3.1.244.1", 1, self.GEN)])
        fields = os.path.join(self.tmp.name, "fields.txt")
        with open(fields, "w") as fh:
            fh.write("3.1.23.1|[1,0,-1,1]\n3.1.244.1|[-6,1,0,1]\n")
        argv = ["find_generators.py", "--degrees", "3", "--method", "thue", "--recheck-generators",
                "--input", fields]
        saved = sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG
        try:
            sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = argv, self.tmp.name, self.log
            with contextlib.redirect_stdout(io.StringIO()) as out, \
                    contextlib.redirect_stderr(io.StringIO()) as err:
                find_generators.main()
        finally:
            sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = saved
        self.assertEqual(read_data(self.path), [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.244.1", 1, self.GEN)])
        self.assertIn("BETTER 3.1.23.1", err.getvalue())
        self.assertIn("2 fields checked", out.getvalue())
        self.assertEqual(self.logged(), [])

    def test_lmfdb_monogenic(self):
        # Only the fields the LMFDB marks 1 and that have no generator here are searched:
        # 3.1.244.1 (a false -1 here, a conflict) and 3.1.31.1 (no row); not 3.1.23.1
        # (already has one), 3.1.44.1 (LMFDB 0) or 3.1.503.1 (LMFDB -1).
        import types
        import find_generators

        lmfdb = [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.31.1", 1, [1, 1, 0, 1]), ("3.1.44.1", 0, [-1, 1, -1, 1]),
                 ("3.1.244.1", 1, [-6, 1, 0, 1]), ("3.1.503.1", -1, [8, 2, -1, 1])]

        class Table:
            def __init__(self):
                self.queries = []

            def rows(self, query):
                self.queries.append(query)
                for label, mono, coeffs in lmfdb:
                    row = {"label": label, "degree": 3, "monogenic": mono, "coeffs": coeffs}
                    if all(row[k] == v for k, v in query.items()):
                        yield row

            def count(self, query):
                return len(list(self.rows(query)))

            def search(self, query, cols):
                return [{c: r[c] for c in cols} for r in self.rows(query)]

        def run(*extra):
            write_data(self.path, [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.244.1", -1, None)])
            table = Table()
            argv = ["find_generators.py", "--degrees", "3", "--method", "pari", "--lmfdb-monogenic"] + list(extra)
            saved = sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG
            sys.modules["lmf"] = types.SimpleNamespace(db=types.SimpleNamespace(nf_fields=table))
            try:
                sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = argv, self.tmp.name, self.log
                with contextlib.redirect_stdout(io.StringIO()) as out, \
                        contextlib.redirect_stderr(io.StringIO()):
                    find_generators.main()
            finally:
                sys.argv, monogenic.DATA_DIR, monogenic.CONFLICT_LOG = saved
                del sys.modules["lmf"]
            self.assertTrue(table.queries)
            self.assertTrue(all(q.get("monogenic") == 1 for q in table.queries))
            return out.getvalue()

        out = run()
        self.assertIn("2 fields checked", out)
        self.assertEqual(read_data(self.path), [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.31.1", 1, [1, 1, 0, 1]),
                                                ("3.1.244.1", 1, self.GEN)])
        self.assertEqual(len(self.logged()), 1)  # the conflict for 3.1.244.1
        out = run("--only-new")
        self.assertIn("1 fields checked", out)
        self.assertEqual(read_data(self.path), [("3.1.23.1", 1, [1, 0, -1, 1]), ("3.1.31.1", 1, [1, 1, 0, 1]),
                                                ("3.1.244.1", -1, None)])


class TestData(unittest.TestCase):
    def test_quadratic(self):
        rows = read_data(os.path.join(DATA_DIR, "degree2.txt"))
        self.assertEqual(rows[:2], [("2.0.3.1", 1, [1, -1, 1]), ("2.0.4.1", 1, [1, 0, 1])])

    def test_header(self):
        import glob
        for path in glob.glob(os.path.join(DATA_DIR, "degree*.txt")):
            with self.subTest(path=os.path.basename(path)), open(path) as fh:
                self.assertEqual(fh.readline(), "label|monogenic|monogenic_coeffs\n")

    def test_old_header(self):
        # Files with the old column name zk_coeffs are still read, and written with the new one.
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "degree2.txt")
            with open(path, "w") as fh:
                fh.write("label|monogenic|zk_coeffs\ntext|smallint|numeric[]\n\n2.0.3.1|1|[1,-1,1]\n")
            rows = read_data(path)
            self.assertEqual(rows, [("2.0.3.1", 1, [1, -1, 1])])
            write_data(path, rows)
            with open(path) as fh:
                self.assertEqual(fh.readline(), "label|monogenic|monogenic_coeffs\n")

    def test_format(self):
        self.assertEqual(format_coeffs([1, 0, -1, 1]), "[1,0,-1,1]")
        self.assertEqual(parse_coeffs("[1,0,-1,1]"), [1, 0, -1, 1])
        self.assertEqual(parse_coeffs("[1, 0, -1, 1]"), [1, 0, -1, 1])


if __name__ == "__main__":
    unittest.main()
