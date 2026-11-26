"""Read number fields from the LMFDB database, and go on when the connection is lost.

A plain `db.nf_fields.search(query)` reads the whole result through one server-side
cursor, so the connection has to stay up for the whole run. After some hours it can
drop (e.g. "SSL error: unexpected eof while reading", then "the connection is lost"),
and the run crashes. search_fields reads the fields a chunk at a time instead, in the
order of |disc| and then of the label, each chunk with a short query of its own: the
next chunk is the fields with |disc| larger than the last one read (and all the fields
with that |disc| together, so none is missed or read twice). Nothing is held open
between chunks, however long the work on a chunk takes. When a query fails because
the connection is lost, it reconnects (db.reset_connection) and runs that query again,
waiting 10 s, 20 s, 40 s, ... (at most 10 min) between tries, so the run goes on where
it was. Other errors (a wrong column, say) are raised at once.

Works with `from lmfdb import db` and with lmfdb-lite (`from lmf import db`).
"""

import sys
import time

CHUNK = 10000    # fields per query
RETRIES = 8      # tries after the first: 10 s + 20 s + ... + 600 s, about 40 min of waiting
WAIT = 10        # seconds before the first retry, doubled each time
MAX_WAIT = 600   # at most this long between two tries


def get_db():
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    return db


def is_connection_error(err):
    """Whether err is a lost or broken connection (psycopg's OperationalError or
    InterfaceError, or a subclass), which is worth reconnecting for."""
    return any(c.__name__ in ("OperationalError", "InterfaceError") for c in type(err).__mro__)


def first_line(err):
    lines = str(err).strip().splitlines()
    return lines[0] if lines else type(err).__name__


def _say(msg):
    print(msg, file=sys.stderr, flush=True)


def with_retries(fetch, db=None, what="reading the LMFDB", say=None, retries=None, wait=None):
    """fetch(), and if that fails because the connection is lost, reconnect and try
    again, up to `retries` more times. fetch should run its query from the start each
    time (e.g. `lambda: list(table.search(...))`)."""
    say = say or _say
    retries = RETRIES if retries is None else retries
    wait = WAIT if wait is None else wait
    for attempt in range(retries + 1):
        try:
            result = fetch()
        except Exception as err:
            if not is_connection_error(err) or attempt == retries:
                raise
            pause = min(wait * 2 ** attempt, MAX_WAIT)
            say("lost the LMFDB connection while %s (%s: %s); reconnecting in %ds (retry %d of %d)"
                % (what, type(err).__name__, first_line(err), pause, attempt + 1, retries))
            time.sleep(pause)
            try:
                (db or get_db()).reset_connection()
            except Exception as err2:  # e.g. the server is still unreachable: try again later
                say("could not reconnect yet (%s: %s)" % (type(err2).__name__, first_line(err2)))
            continue
        if attempt:
            say("reconnected to the LMFDB, going on")
        return result


def count_fields(query, table="nf_fields", db=None, say=None):
    """db.<table>.count(query), reconnecting if the connection is lost."""
    db = db or get_db()
    return with_retries(lambda: getattr(db, table).count(query), db, "counting fields", say)


def label_key(label):
    """Labels n.r.D.i by degree, |disc|, signature, then i, as the LMFDB sorts them."""
    n, r, D, i = (int(x) for x in label.split("."))
    return (n, D, r, i)


def _above(query, lo):
    """query restricted to disc_abs > lo (None if that leaves nothing)."""
    q = dict(query)
    if lo is None:
        return q
    cond = q.get("disc_abs")
    if cond is None:
        q["disc_abs"] = {"$gt": lo}
    elif isinstance(cond, dict):
        cond = dict(cond)
        cond["$gt"] = max(lo, cond["$gt"]) if "$gt" in cond else lo
        q["disc_abs"] = cond
    elif cond <= lo:  # a single discriminant, already read
        return None
    return q


def search_fields(query, cols, table="nf_fields", db=None, chunk=CHUNK, say=None,
                  retries=None, wait=None):
    """The rows of db.<table>.search(query, cols), as dicts (with disc_abs, also when
    it is not in cols), in the order of (disc_abs, label), read `chunk` at a time, and
    reconnecting when the connection is lost (see the top of this file). query may
    restrict disc_abs with $gte/$lte/$gt/$lt or to one value."""
    db = db or get_db()
    cols = list(cols) + ([] if "disc_abs" in cols else ["disc_abs"])

    def fetch(q, limit=None, what="reading fields"):
        if limit is None:
            return with_retries(lambda: list(getattr(db, table).search(q, cols)), db, what, say, retries, wait)
        return with_retries(lambda: list(getattr(db, table).search(q, cols, sort=["disc_abs"], limit=limit)),
                            db, what, say, retries, wait)

    lo = None
    while True:
        q = _above(query, lo)
        if q is None:
            return
        what = "reading fields" + ("" if lo is None else " with |disc| > %s" % lo)
        rows = fetch(q, chunk, what)
        top = None
        if len(rows) >= chunk:
            # The last |disc| in this chunk may go on in the next: read all of it now.
            top = rows[-1]["disc_abs"]
            rows = [r for r in rows if r["disc_abs"] != top]
            rows += fetch(dict(query, disc_abs=top), what="reading the fields with |disc| = %s" % top)
        rows.sort(key=lambda r: label_key(r["label"]))
        yield from rows
        if top is None:
            return
        lo = top
