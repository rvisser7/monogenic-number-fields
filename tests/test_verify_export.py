"""Tests for scripts/verify_export.py. They need Sage (sage -python -m unittest ...) and are
skipped without it."""

import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

try:
    import sage.all  # noqa: F401
    from verify_export import check, label_disc
except ImportError:
    sage = None

# The two cyclic cubic fields of conductor 63, both of discriminant 3969 (nflist("C3", 3969)).
# x^3 - 21x - 35 has discriminant 3969, so its root generates O_K; x^3 - 21x - 28 has
# discriminant 4 * 3969.
C63_1 = [-28, -21, 0, 1]
C63_2 = [-35, -21, 0, 1]


@unittest.skipIf(sage is None, "needs Sage")
class TestVerifyExport(unittest.TestCase):
    def test_label_disc(self):
        self.assertEqual(label_disc("2.0.3.1"), -3)
        self.assertEqual(label_disc("4.2.1600.1"), -1600)
        self.assertEqual(label_disc("3.3.49.1"), 49)

    def test_good(self):
        self.assertIsNone(check("1.1.1.1", [0, 1], ([0, 1], 1)))
        self.assertIsNone(check("2.0.3.1", [1, -1, 1], ([1, -1, 1], -3)))
        # another generator of the same field: alpha = 1 - zeta_3
        self.assertIsNone(check("2.0.3.1", [3, -3, 1], ([1, -1, 1], -3)))
        self.assertIsNone(check("3.3.49.1", [1, -2, -1, 1], ([1, -2, -1, 1], 49)))
        self.assertIsNone(check("3.3.49.1", [1, -2, -1, 1], None))

    def test_bad(self):
        cases = [
            ("2.0.3.1", [1, -1, 2], ([1, -1, 1], -3), "monic"),
            ("2.0.3.1", [1, -1, 0, 1], ([1, -1, 1], -3), "degree"),
            ("2.2.4.1", [-1, 0, 1], ([-1, 0, 1], 4), "reducible"),
            ("2.0.4.1", [5, 0, 1], ([1, 0, 1], -4), "disc(f)"),       # index 2
            ("2.0.4.1", [1, -1, 1], ([1, 0, 1], -4), "disc(f)"),      # another field
            ("2.0.3.1", [1, -1, 1], ([1, -1, 1], -4), "disc(K)"),     # LMFDB disagrees
            ("3.3.3969.1", C63_2, (C63_1, 3969), "not isomorphic"),   # the other field
            ("3.3.3969.1", C63_1, (C63_1, 3969), "disc(f)"),          # index 2
        ]
        for label, coeffs, field, reason in cases:
            with self.subTest(label=label, coeffs=coeffs):
                self.assertIn(reason, check(label, coeffs, field) or "")

    def test_same_disc_fields(self):
        self.assertIsNone(check("3.3.3969.2", C63_2, (C63_2, 3969)))
        self.assertIsNone(check("3.3.3969.2", [35, -21, 0, 1], (C63_2, 3969)))  # -alpha

    def test_command_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            export = os.path.join(tmp, "export.txt")
            with open(export, "w") as fh:
                fh.write("label|monogenic_coeffs\ntext|numeric[]\n\n"
                         "2.0.3.1|[1,-1,1]\n2.0.4.1|[5,0,1]\n")
            fields = os.path.join(tmp, "fields.txt")
            with open(fields, "w") as fh:
                fh.write("2.0.3.1|[1,-1,1]\n2.0.4.1|[1,0,1]\n")
            run = subprocess.run([sys.executable, os.path.join(SCRIPTS, "verify_export.py"),
                                  export, "--fields", fields], capture_output=True, text=True)
            self.assertEqual(run.returncode, 1)
            self.assertEqual([line.split("|")[0] for line in run.stdout.split("\n") if line],
                             ["2.0.4.1"])
            self.assertIn("2 rows checked, 1 failed", run.stderr)


if __name__ == "__main__":
    unittest.main()
