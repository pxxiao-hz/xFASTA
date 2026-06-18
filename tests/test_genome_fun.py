import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genome_fun import build_genomescope_commands, genome_split_chr_by_gap  # noqa: E402


class GenomeFunTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.workdir = Path(self.tmpdir.name)
        self.previous_cwd = Path.cwd()
        os.chdir(self.workdir)

    def tearDown(self):
        os.chdir(self.previous_cwd)

    def test_split_chr_by_gap_splits_on_any_n_by_default(self):
        fasta = self.workdir / "genome.fa"
        fasta.write_text(">chr1\nAAAANNCCNNTT\n>chr2\nNNGG\n")

        genome_split_chr_by_gap(str(fasta))

        self.assertEqual(
            (self.workdir / "genome.contig.fa").read_text(),
            ">chr1_contig1\nAAAA\n"
            ">chr1_contig2\nCC\n"
            ">chr1_contig3\nTT\n"
            ">chr2_contig1\nGG\n",
        )

    def test_split_chr_by_gap_respects_minimum_gap_size(self):
        fasta = self.workdir / "genome.fa"
        fasta.write_text(">chr1\nAAAANCCNNNTT\n")

        genome_split_chr_by_gap(str(fasta), 2)

        self.assertEqual(
            (self.workdir / "genome.contig.fa").read_text(),
            ">chr1_contig1\nAAAANCC\n"
            ">chr1_contig2\nTT\n",
        )

    def test_build_genomescope_commands_uses_configured_scripts(self):
        cmd_v1, cmd_v2 = build_genomescope_commands(
            21,
            150,
            4,
            "/opt/genomescope1/genomescope.R",
            "/opt/genomescope2/genomescope.R",
        )

        self.assertEqual(
            cmd_v1,
            "Rscript /opt/genomescope1/genomescope.R reads.21.histo 21 150 21_genomescope1_output",
        )
        self.assertEqual(
            cmd_v2,
            "Rscript /opt/genomescope2/genomescope.R -i reads.21.histo -o 21_genomescope2_output -k 21 -p 4",
        )


if __name__ == "__main__":
    unittest.main()
