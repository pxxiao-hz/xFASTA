import gzip
import contextlib
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from reads_fun import (  # noqa: E402
    build_chopper_filter_command,
    func_get_reads_phreads,
    func_reads_ID_simplified,
    func_reads_deduplication_based_ID,
    func_reads_extract_based_ID,
    func_reads_filter_ONT_reads,
)


FASTQ_TEXT = """@read1 old description
ACGT
+
IIII
@read1 duplicate description
TGCA
+
!!!!
@read2 second description
AAAA
+
####
"""


class ReadsFunTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.workdir = Path(self.tmpdir.name)
        self.previous_cwd = Path.cwd()
        os.chdir(self.workdir)

    def tearDown(self):
        os.chdir(self.previous_cwd)

    def write_fastq(self, name="reads.fastq.gz"):
        path = self.workdir / name
        if name.endswith(".gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(FASTQ_TEXT)
        else:
            path.write_text(FASTQ_TEXT)
        return path

    def read_gzip_text(self, path):
        with gzip.open(path, "rt") as handle:
            return handle.read()

    def test_deduplication_keeps_first_record_for_each_read_id(self):
        input_fastq = self.write_fastq()

        func_reads_deduplication_based_ID(str(input_fastq))

        output = self.read_gzip_text(self.workdir / "reads.redu.fq.gz")
        self.assertEqual(output.count("@read1"), 1)
        self.assertIn("@read2", output)
        self.assertIn("ACGT", output)
        self.assertNotIn("TGCA", output)

    def test_id_simplified_streams_each_record_once(self):
        input_fastq = self.write_fastq()

        func_reads_ID_simplified(str(input_fastq))

        output = self.read_gzip_text(self.workdir / "reads.simplified.fastq.gz")
        self.assertEqual(output.count("@read1"), 2)
        self.assertEqual(output.count("@read2"), 1)
        self.assertNotIn("old description", output)
        self.assertNotIn("duplicate description", output)

    def test_extract_supports_plain_fastq_with_full_input_path(self):
        input_fastq = self.write_fastq("nested.fastq")

        func_reads_extract_based_ID(str(input_fastq), "read2")

        output = self.read_gzip_text(self.workdir / "nested.extract.fastq.gz")
        self.assertIn("@read2", output)
        self.assertIn("AAAA", output)
        self.assertNotIn("@read1", output)

    def test_phreads_uses_input_path_and_supports_gzip(self):
        input_fastq = self.write_fastq()

        with contextlib.redirect_stdout(io.StringIO()):
            func_get_reads_phreads(str(input_fastq))

        output = (self.workdir / "reads.qv.txt").read_text().splitlines()
        self.assertEqual(len(output), 3)
        self.assertEqual(float(output[0]), 40.0)

    def test_build_chopper_filter_command_returns_argument_list(self):
        command, output = build_chopper_filter_command(
            "reads dir/sample.fastq.gz",
            7,
            100000,
            "/opt/chopper bin/chopper",
        )

        self.assertEqual(output, "sample.q7.l100000.fastq.gz")
        self.assertEqual(
            command,
            ["/opt/chopper bin/chopper", "-q", "7", "-l", "100000"],
        )

    def test_filter_ont_uses_pipeline_without_shell(self):
        input_fastq = self.write_fastq()
        fake_chopper = self.workdir / "fake-chopper"
        fake_chopper.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "sys.stdout.buffer.write(sys.stdin.buffer.read())\n"
        )
        fake_chopper.chmod(0o755)

        func_reads_filter_ONT_reads(str(input_fastq), 7, 3, str(fake_chopper))

        output = self.read_gzip_text(self.workdir / "reads.q7.l3.fastq.gz")
        self.assertEqual(output, FASTQ_TEXT)


if __name__ == "__main__":
    unittest.main()
