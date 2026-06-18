import contextlib
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gap_fun import func_get_gap_location  # noqa: E402


class GapFunTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.workdir = Path(self.tmpdir.name)
        self.previous_cwd = Path.cwd()
        os.chdir(self.workdir)

    def tearDown(self):
        os.chdir(self.previous_cwd)

    def test_gap_location_writes_gff_and_half_open_bed_coordinates(self):
        fasta = self.workdir / "genome.fa"
        fasta.write_text(">chr1\nAANNNCCnN\n>chr2\nNNAA\n")

        with contextlib.redirect_stdout(io.StringIO()):
            func_get_gap_location(str(fasta))

        self.assertEqual(
            (self.workdir / "genome.gaps.gff").read_text(),
            "chr1\t.\tgap\t3\t5\t.\t.\t.\tName=gap1;size=3\n"
            "chr1\t.\tgap\t8\t9\t.\t.\t.\tName=gap2;size=2\n"
            "chr2\t.\tgap\t1\t2\t.\t.\t.\tName=gap3;size=2\n",
        )
        self.assertEqual(
            (self.workdir / "genome.gaps.bed").read_text(),
            "chr1\t2\t5\n"
            "chr1\t7\t9\n"
            "chr2\t0\t2\n",
        )


if __name__ == "__main__":
    unittest.main()
