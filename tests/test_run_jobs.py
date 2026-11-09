import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

from run_jobs import log_name, read_jobs, run_jobs  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestRunJobs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.scripts = os.path.join(self.dir, "scripts")
        self.logs = os.path.join(self.dir, "logs")
        os.makedirs(self.scripts)
        # A stand-in for find_generators.py: prints its arguments, fails if asked to.
        with open(os.path.join(self.scripts, "job.py"), "w") as fh:
            fh.write("import sys, time\nprint('args', sys.argv[1:])\ntime.sleep(0.3)\n"
                     "sys.exit(3 if '--fail' in sys.argv else 0)\n")

    def tearDown(self):
        self.tmp.cleanup()

    def write_jobs(self, text):
        path = os.path.join(self.dir, "jobs.txt")
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_read_jobs(self):
        path = self.write_jobs("# a comment\n\nfind_generators.py --degrees 5 -B 4  # five\n"
                               "  find_generators.py --degrees 6 -B 3\n")
        self.assertEqual(read_jobs(path), [["find_generators.py", "--degrees", "5", "-B", "4"],
                                           ["find_generators.py", "--degrees", "6", "-B", "3"]])

    def test_log_name(self):
        self.assertEqual(log_name(3, ["find_generators.py", "--degrees", "5", "-B", "4"]),
                         "03-find_generators-d5.log")
        self.assertEqual(log_name(12, ["find_generators.py", "--degrees=5-8"]), "12-find_generators-d5-8.log")
        self.assertEqual(log_name(1, ["check_data.py"]), "01-check_data.log")

    def test_run(self):
        jobs = [["job.py", "--degrees", "5"], ["job.py", "--degrees", "6", "--fail"], ["job.py", "--degrees", "7"]]
        out = io.StringIO()
        codes = run_jobs(jobs, 2, self.logs, self.scripts, out)
        self.assertEqual(codes, [0, 3, 0])
        self.assertEqual(sorted(os.listdir(self.logs)),
                         ["01-job-d5.log", "02-job-d6.log", "03-job-d7.log"])
        with open(os.path.join(self.logs, "02-job-d6.log")) as fh:
            self.assertEqual(fh.read(), "args ['--degrees', '6', '--fail']\n")
        text = out.getvalue()
        self.assertEqual(text.count("started"), 3)
        self.assertIn("FAILED (exit 3) 2/3", text)
        self.assertIn("all 3 jobs finished, 1 failed: 2", text)

    def test_missing_script(self):
        with self.assertRaises(SystemExit):
            run_jobs([["nope.py"]], 1, self.logs, self.scripts, io.StringIO())

    def test_example_jobs_file(self):
        # Every job in jobs.txt names a script that exists.
        for job in read_jobs(os.path.join(REPO, "jobs.txt")):
            self.assertTrue(os.path.isfile(os.path.join(REPO, "scripts", job[0])), job[0])


if __name__ == "__main__":
    unittest.main()
