import os
import random
import sys
import tempfile
import unittest

TESTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(TESTS), "scripts"))
sys.path.insert(0, TESTS)

from fake_lmfdb import FakeDB, OperationalError  # noqa: E402
from lmfdb_fields import count_fields, label_key, search_fields, with_retries  # noqa: E402


def fields(seed=1, n=500):
    """Degree 3 and 4 fields, many sharing |disc| (up to 30 for one |disc|, more than a chunk)."""
    rng = random.Random(seed)
    out = []
    for deg in (3, 4):
        k = {}
        for _ in range(n):
            D = rng.choice(range(23, 23 + n // 5))
            r = rng.choice([1, 3]) if deg == 3 else rng.choice([0, 2, 4])
            k[r, D] = k.get((r, D), 0) + 1
            out.append({"label": "%d.%d.%d.%d" % (deg, r, D, k[r, D]), "coeffs": [D, r, 1]})
    out += [{"label": "3.1.10000.%d" % i, "coeffs": [i]} for i in range(1, 31)]
    rng.shuffle(out)
    return out


class TestSearchFields(unittest.TestCase):
    def expected(self, rows, query):
        db = FakeDB(nf_fields=rows)
        return sorted((r["label"] for r in rows if db.nf_fields._match(dict(r, degree=int(r["label"][0]),
                                                                             disc_abs=int(r["label"].split(".")[2])),
                                                                       query)), key=label_key)

    def read(self, db, query, chunk, said=None):
        return [r["label"] for r in search_fields(query, ["label", "coeffs"], db=db, chunk=chunk, wait=0,
                                                  say=(said.append if said is not None else lambda m: None))]

    def test_order_and_chunks(self):
        rows = fields()
        for query in ({"degree": 3}, {"degree": 4, "disc_abs": {"$gte": 40, "$lte": 90}},
                      {"degree": 3, "disc_abs": 10000}, {"degree": 3, "disc_abs": {"$gt": 50}}):
            want = self.expected(rows, query)
            self.assertTrue(want)
            for chunk in (1, 2, 7, 100, 10**4):
                with self.subTest(query=query, chunk=chunk):
                    got = self.read(FakeDB(nf_fields=rows), query, chunk)
                    self.assertEqual(got, want)

    def test_lost_connection(self):
        # The connection drops in the middle of chunks, during the read of all the fields
        # with one |disc|, and twice in a row: every field comes once, in order.
        rows = fields(2)
        want = self.expected(rows, {"degree": 3})
        for drops in ([5], [5, 6], [1], [40, 41, 42, 300], list(range(10, 400, 37))):
            with self.subTest(drops=drops):
                db = FakeDB(drop_at=drops, nf_fields=rows)
                said = []
                self.assertEqual(self.read(db, {"degree": 3}, 7, said), want)
                self.assertEqual(db.resets, len(drops))
                self.assertTrue(any("lost the LMFDB connection" in m for m in said))
                self.assertTrue(any("reconnected" in m for m in said))

    def test_lazy(self):
        # Nothing is read before it is needed, and a chunk is read in one go.
        db = FakeDB(nf_fields=fields())
        it = search_fields({"degree": 3}, ["label"], db=db, chunk=10)
        self.assertEqual(db.uses, 0)
        next(it)
        self.assertLessEqual(db.uses, 40)  # the chunk, and the fields with its last |disc|

    def test_count(self):
        db = FakeDB(drop_at=[1], nf_fields=fields())
        self.assertEqual(count_fields({"degree": 3}, db=db, say=lambda m: None),
                         sum(1 for r in fields() if r["label"][0] == "3"))
        self.assertEqual(db.resets, 1)

    def test_gives_up(self):
        db = FakeDB(nf_fields=fields())
        db.lost = True

        def reset():  # the server stays down
            db.resets += 1
        db.reset_connection = reset
        with self.assertRaises(OperationalError):
            list(search_fields({"degree": 3}, ["label"], db=db, wait=0, retries=3, say=lambda m: None))
        self.assertEqual(db.resets, 3)

    def test_other_errors(self):
        # Not a connection problem: raised at once, no retry.
        calls = []

        def fetch():
            calls.append(1)
            raise KeyError("no such column")
        with self.assertRaises(KeyError):
            with_retries(fetch, FakeDB(), wait=0, say=lambda m: None)
        self.assertEqual(len(calls), 1)


class TestFindGenerators(unittest.TestCase):
    def test_lost_connection(self):
        # find_generators.py reads every field once, in order, though the connection
        # drops again and again.
        try:
            import cypari2
        except ImportError:
            self.skipTest("needs cypari2")
        import contextlib
        import io
        import types
        from unittest import mock
        import find_generators
        pari = cypari2.Pari()
        seen, rows = {}, []
        for a in range(-6, 7):
            for b in range(1, 8):
                f = pari("x^3 + %d*x + %d" % (a, b))
                if not f.polisirreducible():
                    continue
                f = f.polredabs()
                if str(f) in seen:
                    continue
                seen[str(f)] = 1
                D = int(pari.nfdisc(f))
                r = int(pari.polsturm(f))
                k = sum(1 for lab in (x["label"] for x in rows) if lab.startswith("3.%d.%d." % (r, abs(D)))) + 1
                rows.append({"label": "3.%d.%d.%d" % (r, abs(D), k), "coeffs": [int(c) for c in reversed(f.Vec())]})

        def run(drops):
            db = FakeDB(drop_at=drops, nf_fields=rows)
            argv = ["find_generators.py", "--degrees", "3", "--method", "polredabs", "--print", "--chunk", "4"]
            saved = sys.argv
            sys.modules["lmf"] = types.SimpleNamespace(db=db)
            try:
                sys.argv = argv
                with contextlib.redirect_stdout(io.StringIO()) as out, \
                        contextlib.redirect_stderr(io.StringIO()) as err, mock.patch("lmfdb_fields.WAIT", 0), \
                        tempfile.TemporaryDirectory() as tmp, mock.patch("monogenic.DATA_DIR", tmp):
                    find_generators.main()
            finally:
                sys.argv = saved
                del sys.modules["lmf"]
            return out.getvalue(), err.getvalue(), db.resets

        out, err, resets = run([])
        self.assertEqual(resets, 0)
        self.assertGreater(out.count("|1|"), 20)
        self.assertIn("%d fields checked" % len(rows), err)
        out2, err2, resets = run([1, 7, 8, 30, 31, 32, 50])
        self.assertEqual(resets, 7)
        self.assertEqual(out2, out)
        self.assertIn("%d fields checked" % len(rows), err2)
        self.assertIn("lost the LMFDB connection", err2)


if __name__ == "__main__":
    unittest.main()
