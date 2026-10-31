import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

from monogenic import DATA_DIR, find_generator, format_coeffs, parse_coeffs, label_disc, pari, read_data  # noqa: E402

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

    def test_dedekind_non_monogenic(self):
        # Dedekind's cubic x^3 - x^2 - 2x - 8 is not monogenic (2 splits completely).
        self.assertIsNone(find_generator([-8, -2, -1, 1]))


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
