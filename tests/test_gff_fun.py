import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gff_fun import func_gff_NewGFF  # noqa: E402


class GffFunTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.workdir = Path(self.tmpdir.name)
        self.previous_cwd = Path.cwd()
        os.chdir(self.workdir)

    def tearDown(self):
        os.chdir(self.previous_cwd)

    def test_new_gff_generates_missing_gene_and_exon_without_shell(self):
        gff = self.workdir / "input.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1\n"
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1\n"
            "chr1\tsrc\tmisc\t1\t30\t.\t+\t.\tID=ignored\n"
        )

        func_gff_NewGFF(str(gff))

        self.assertEqual(
            (self.workdir / "input.new.gff").read_text(),
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "chr1\tsrc\texon\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1\n"
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1\n",
        )


if __name__ == "__main__":
    unittest.main()
