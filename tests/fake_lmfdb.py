"""A stand-in for the LMFDB database (psycodict), enough for scripts/lmfdb_fields.py:
tables with count and search (query with $ne/$gt/$gte/$lt/$lte/$in, sort, limit), and
a database whose connection can be made to drop."""


class OperationalError(Exception):
    """Named like psycopg's, which is what lmfdb_fields reconnects for."""


class FakeTable:
    def __init__(self, rows, db=None):
        self.rows = []
        for r in rows:
            r = dict(r)
            n, _, D, _ = (int(x) for x in r["label"].split("."))
            r.setdefault("degree", n)
            r.setdefault("disc_abs", D)
            self.rows.append(r)
        self.db = db
        self.queries = []

    def _match(self, row, query):
        for col, cond in query.items():
            val = row.get(col)
            if isinstance(cond, dict):
                for op, x in cond.items():
                    if op == "$in":
                        if val not in x:
                            return False
                        continue
                    if val is None:  # SQL: comparisons with NULL are not true
                        return False
                    if not {"$ne": val != x, "$gt": val > x, "$gte": val >= x,
                            "$lt": val < x, "$lte": val <= x}[op]:
                        return False
            elif val != cond:
                return False
        return True

    def count(self, query):
        if self.db:
            self.db.use()
        return sum(1 for r in self.rows if self._match(r, query))

    def search(self, query, cols, sort=None, limit=None):
        self.queries.append(query)
        rows = [r for r in self.rows if self._match(r, query)]
        if sort:
            rows.sort(key=lambda r: tuple(r[c] for c in sort))
        if limit is not None:
            rows = rows[:limit]
        for r in rows:
            if self.db:
                self.db.use()
            yield {c: r[c] for c in cols if r.get(c) is not None}


class FakeDB:
    """Tables as attributes. drop_at: the numbers of the uses (rows read or counts) at
    which the connection is lost; it stays lost until reset_connection()."""

    def __init__(self, drop_at=(), **tables):
        self.uses = 0
        self.drop_at = set(drop_at)
        self.lost = False
        self.resets = 0
        for name, rows in tables.items():
            setattr(self, name, FakeTable(rows, self))

    def use(self):
        if self.lost:
            raise OperationalError("the connection is lost")
        self.uses += 1
        if self.uses in self.drop_at:
            self.lost = True
            raise OperationalError("consuming input failed: SSL error: unexpected eof while reading")

    def reset_connection(self):
        self.resets += 1
        self.lost = False
