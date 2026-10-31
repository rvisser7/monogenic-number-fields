"""Progress reporting on stderr for the long-running scripts."""

import sys
import time


class Progress:
    """Progress on stderr: one line rewritten in place on a terminal, otherwise
    a line every `every` seconds (so log files stay readable)."""

    def __init__(self, every=30):
        self.tty = sys.stderr.isatty()
        self.every = 0.2 if self.tty else every
        self.last = 0

    def start(self, what):
        self.what = what
        self.t0 = time.time()
        self.last = 0

    def __call__(self, done, total, extra=""):
        now = time.time()
        if now - self.last < self.every and done < total:
            return
        self.last = now
        elapsed = now - self.t0
        eta = elapsed / done * (total - done) if done else 0
        line = (f"{self.what}: {done}/{total} ({100 * done // total}%){extra}, "
                f"{fmt_time(elapsed)} elapsed, ~{fmt_time(eta)} left")
        if self.tty:
            print("\r\033[K" + line, end="\n" if done == total else "", file=sys.stderr, flush=True)
        else:
            print(line, file=sys.stderr, flush=True)


def fmt_time(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m{s:02d}s" if h else f"{m}m{s:02d}s" if m else f"{s}s"
