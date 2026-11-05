import contextlib
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import monogenic  # noqa: E402
from monogenic import (DATA_DIR, cubic_generator, find_generator, format_coeffs, label_disc,  # noqa: E402
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
        self.assertEqual(find_generator(ISSUE_EXAMPLES["4.0.1040.2"]), [10, -10, 3, -2, 1])

    def test_root_generates(self):
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
        other = [10, -6, 1, 1]
        self.assertEqual(self.merge((1, self.GEN), (1, other))[:2], ((1, self.GEN), (0, 0, 1, 0)))
        self.assertEqual(self.merge((1, self.GEN), (1, other), True)[:2], ((1, other), (0, 1, 0, 0)))
        self.assertEqual(self.logged(), [])

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


class TestData(unittest.TestCase):
    def test_quadratic(self):
        rows = read_data(os.path.join(DATA_DIR, "degree2.txt"))
        self.assertEqual(rows[:2], [("2.0.3.1", 1, [1, -1, 1]), ("2.0.4.1", 1, [1, 0, 1])])

    def test_format(self):
        self.assertEqual(format_coeffs([1, 0, -1, 1]), "[1,0,-1,1]")
        self.assertEqual(parse_coeffs("[1,0,-1,1]"), [1, 0, -1, 1])
        self.assertEqual(parse_coeffs("[1, 0, -1, 1]"), [1, 0, -1, 1])


if __name__ == "__main__":
    unittest.main()
