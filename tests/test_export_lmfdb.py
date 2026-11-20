import io
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

from export_lmfdb import export  # noqa: E402

HEADER = "label|monogenic|monogenic_coeffs\ntext|smallint|numeric[]\n\n"


class TestExportLmfdb(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, rows, header=HEADER):
        path = os.path.join(self.dir, name)
        with open(path, "w") as fh:
            fh.write(header + "".join(row + "\n" for row in rows))
        return path

    def test_export(self):
        paths = [
            self.write("degree10.txt", ["10.0.1000.1|1|[1,0,0,0,0,0,0,0,0,0,1]"]),
            self.write("degree2.txt", ["2.2.5.1|1|[-1,-1,1]", "2.0.3.1|1|[1,-1,1]",
                                       "2.2.40.1|-1|\\N", "2.0.4.1|1|[1,0,1]",
                                       "2.2.12.1|0|\\N"]),
            self.write("degree1.txt", ["1.1.1.1|1|[0,1]"]),
            self.write("degree3.txt", ["3.1.23.1|1|[-1,-1,0,1]"],
                       header=HEADER.replace("monogenic_coeffs", "zk_coeffs")),
        ]
        out = io.StringIO()
        self.assertEqual(export(paths, out), 6)
        self.assertEqual(out.getvalue(),
                         "label|monogenic_coeffs\ntext|numeric[]\n\n"
                         "1.1.1.1|[0,1]\n2.0.3.1|[1,-1,1]\n2.0.4.1|[1,0,1]\n2.2.5.1|[-1,-1,1]\n"
                         "3.1.23.1|[-1,-1,0,1]\n10.0.1000.1|[1,0,0,0,0,0,0,0,0,0,1]\n")

    def test_bad_rows(self):
        for row in ["2.0.3.1|1|\\N", "2.0.3.1|1|[1,-1]", "2.0.3.1|1|[1, -1, 1]"]:
            path = self.write("degree2.txt", [row])
            with self.assertRaises(ValueError):
                export([path], io.StringIO())
        path = self.write("degree2.txt", ["2.0.3.1|1|[1,-1,1]"], header="label|x\n\n\n")
        with self.assertRaises(ValueError):
            export([path], io.StringIO())

    def test_command_line(self):
        self.write("degree2.txt", ["2.0.3.1|1|[1,-1,1]", "2.0.4.1|-1|\\N"])
        self.write("degree3.txt", ["3.1.23.1|1|[-1,-1,0,1]"])
        out = os.path.join(self.dir, "out.txt")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "export_lmfdb.py"),
                        "--data-dir", self.dir, "--degrees", "2", "--output", out],
                       check=True, capture_output=True)
        with open(out) as fh:
            self.assertEqual(fh.read(), "label|monogenic_coeffs\ntext|numeric[]\n\n"
                                        "2.0.3.1|[1,-1,1]\n")


if __name__ == "__main__":
    unittest.main()
