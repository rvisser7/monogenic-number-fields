import os
import subprocess
import sys
import tempfile
import types
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

import compare_lmfdb  # noqa: E402
from compare_lmfdb import DBSource, FileSource, Log, compare  # noqa: E402
import contextlib  # noqa: E402
import io  # noqa: E402

# label: (LMFDB monogenic, LMFDB coeffs), as in nf_fields (monogenic may be null)
LMFDB = {
    "3.1.23.1": (1, [-1, -1, 0, 1]),     # agree 1
    "3.1.31.1": (0, [-1, 1, 0, 1]),      # only ours 1
    "3.1.44.1": (None, [-1, 1, -1, 1]),  # only ours 1 (null)
    "3.1.76.1": (-1, [-2, -2, 0, 1]),    # conflict: ours 1
    "3.1.104.1": (1, [-2, 1, -1, 1]),    # only LMFDB 1
    "3.1.107.1": (-1, [-3, 1, -1, 1]),   # only LMFDB -1
    "3.3.49.1": (2, [1, -2, -1, 1]),     # unexpected value
    "3.1.503.1": (-1, [-4, -1, -1, 1]),  # agree -1
    "3.1.588.1": (1, [-8, 0, 0, 1]),     # conflict: ours -1
    "3.3.81.1": (0, [-1, -3, 0, 1]),     # only ours -1
    "4.0.117.1": (1, [3, -1, 0, 0, 1]),  # other degree
}

OURS = {
    "3.1.23.1": (1, [-1, -1, 0, 1]),
    "3.1.31.1": (1, [-1, 1, 0, 1]),
    "3.1.44.1": (1, [-1, 1, -1, 1]),
    "3.1.76.1": (1, [-2, -2, 0, 1]),
    "3.1.503.1": (-1, None),
    "3.1.588.1": (-1, None),
    "3.3.81.1": (-1, None),
    "3.1.200.1": (0, None),              # explicit 0 row, LMFDB has nothing
}


class FakeTable:
    """Enough of a psycodict table for DBSource."""

    def __init__(self, rows):
        self.rows = [{"label": label, "degree": int(label.split(".")[0]), "disc_abs": int(label.split(".")[2]),
                      "monogenic": mono, "coeffs": coeffs}
                     for label, (mono, coeffs) in rows.items()]

    def _match(self, row, query):
        for col, cond in query.items():
            if isinstance(cond, dict):
                for op, val in cond.items():
                    if row[col] is None:  # SQL: NULL != 0 is not true
                        return False
                    if not {"$ne": row[col] != val, "$gte": row[col] >= val, "$lte": row[col] <= val}[op]:
                        return False
            elif row[col] != cond:
                return False
        return True

    def search(self, query, cols):
        return [{c: r[c] for c in cols} for r in self.rows if self._match(r, query)]

    def count(self, query):
        return sum(1 for r in self.rows if self._match(r, query))


class TestCompare(unittest.TestCase):
    def check(self, result):
        self.assertEqual(result["agree_yes"], ["3.1.23.1"])
        self.assertEqual(result["agree_no"], ["3.1.503.1"])
        self.assertEqual(result["conflict"], ["3.1.76.1", "3.1.588.1"])
        self.assertEqual(result["ours_only_yes"], ["3.1.31.1", "3.1.44.1"])
        self.assertEqual(result["ours_only_no"], ["3.3.81.1"])
        self.assertEqual(result["lmfdb_only_yes"], ["3.1.104.1"])
        self.assertEqual(result["lmfdb_only_no"], ["3.1.107.1"])
        self.assertEqual(result["unexpected"], ["3.3.49.1"])

    def test_db_source(self):
        sys.modules["lmf"] = types.SimpleNamespace(db=types.SimpleNamespace(nf_fields=FakeTable(LMFDB)))
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                known, counts = DBSource(with_coeffs=True).degree(3, Log(every=0))
        finally:
            del sys.modules["lmf"]
        # progress: the counts, then a line per row read (every=0), then the total
        lines = out.getvalue().splitlines()
        self.assertIn("degree 3: 10 LMFDB fields, 2 with monogenic = 0", lines[1])
        self.assertIn("1 LMFDB fields with monogenic != 0 read so far", lines[2])
        self.assertIn("degree 3: 7 LMFDB fields with monogenic != 0 read", lines[-1])
        self.assertEqual(set(known), {label for label, (m, _) in LMFDB.items() if m not in (0, None) and label[0] == "3"})
        self.assertEqual(known["3.1.104.1"], (1, [-2, 1, -1, 1]))
        self.assertEqual((counts[1], counts[-1], counts[0], counts[None], counts["other"], counts["total"]),
                         (3, 3, 2, 1, 1, 10))
        self.check(compare(OURS, known))

    def test_file_source(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write("label|monogenic|coeffs\ntext|smallint|numeric[]\n\n")
            for label, (mono, coeffs) in LMFDB.items():
                fh.write("%s|%s|%s\n" % (label, r"\N" if mono is None else mono,
                                         "[" + ",".join(map(str, coeffs)) + "]"))
        try:
            source = FileSource(fh.name)
        finally:
            os.unlink(fh.name)
        self.assertEqual(source.degrees(), [3, 4])
        known, counts = source.degree(3, Log(quiet=True))
        self.assertEqual((counts[1], counts[-1], counts[0], counts[None], counts["other"], counts["total"]),
                         (3, 3, 2, 1, 1, 10))
        self.check(compare(OURS, known))

    def test_disc_range(self):
        from compare_lmfdb import DiscRange, parse_int
        self.assertEqual(parse_int("10^4"), 10000)
        self.assertEqual(parse_int("1e4"), 10000)
        disc = DiscRange(None, 100)
        self.assertEqual(disc.query(), {"disc_abs": {"$lte": 100}})
        self.assertEqual(DiscRange().query(), {})
        self.assertIn("3.1.76.1", disc)
        self.assertNotIn("3.1.104.1", disc)
        self.assertIn("3.1.104.1", DiscRange(104, 104))
        # The database query gets the disc_abs condition; FakeTable checks it.
        sys.modules["lmf"] = types.SimpleNamespace(db=types.SimpleNamespace(nf_fields=FakeTable(LMFDB)))
        try:
            known, counts = DBSource(False, DiscRange(50, 200)).degree(3, Log(quiet=True))
        finally:
            del sys.modules["lmf"]
        self.assertEqual(sorted(known), ["3.1.104.1", "3.1.107.1", "3.1.76.1"])
        self.assertEqual(counts["total"], 4)  # and 3.3.81.1, with monogenic = 0

    def test_script_on_data(self):
        # Against a copy of our own degree 3 rows the script finds no conflicts, and
        # flipping one row gives a conflict and exit status 1.
        rows = compare_lmfdb.read_data(compare_lmfdb.data_path(3))[:200]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "lmfdb.txt")
            with open(path, "w") as fh:
                fh.write("".join("%s|%d\n" % (label, mono) for label, mono, _ in rows))
            run = [sys.executable, os.path.join(SCRIPTS, "compare_lmfdb.py"), "--lmfdb", path, "--degrees", "3"]
            proc = subprocess.run(run, capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("0 conflicts", proc.stdout)
            self.assertIn("200 LMFDB fields; 200 agree, 0 conflicts", proc.stdout)  # progress line
            label, mono, _ = rows[0]
            with open(path, "a") as fh:
                fh.write("%s|%d\n" % (label, -mono))  # the later line wins
            proc = subprocess.run(run + ["--write", tmp], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 1, proc.stdout + proc.stderr)
            self.assertIn("CONFLICT " + label, proc.stdout)
            with open(os.path.join(tmp, "conflicts.txt")) as fh:
                self.assertTrue(fh.read().startswith("%s|%d|%d|" % (label, mono, -mono)))


if __name__ == "__main__":
    unittest.main()
