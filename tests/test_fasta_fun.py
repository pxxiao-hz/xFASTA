import contextlib
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fasta_fun import (  # noqa: E402
    extract_based_coords,
    extract_based_keywords,
    extract_based_length,
    fasta_extract_based_length,
    func_fasta_compare,
    func_get_fasta_length,
    genome_change_chr_name,
    genome_karyotype,
    genome_reverse_some_chr,
)
from util import split_fasta_based_bp, split_fasta_based_number  # noqa: E402


FASTA_TEXT = """>chr1 description
ACGTACGTAC
GT
>chr2
NNNNAA
>scaffold_3
TTTT
"""


class FastaFunTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.workdir = Path(self.tmpdir.name)
        self.previous_cwd = Path.cwd()
        os.chdir(self.workdir)
        self.fasta = self.workdir / "genome.fa"
        self.fasta.write_text(FASTA_TEXT)

    def tearDown(self):
        os.chdir(self.previous_cwd)

    def test_split_fasta_based_bp_streams_each_record(self):
        with contextlib.redirect_stdout(io.StringIO()):
            split_fasta_based_bp(str(self.fasta), 5)

        self.assertEqual((self.workdir / "chr1 description_1.fa").read_text(), ">chr1 description_1\nACGTA\n")
        self.assertEqual((self.workdir / "chr1 description_3.fa").read_text(), ">chr1 description_3\nGT\n")
        self.assertEqual((self.workdir / "chr2_2.fa").read_text(), ">chr2_2\nA\n")

    def test_split_fasta_based_number_writes_requested_records_per_file(self):
        split_fasta_based_number(str(self.fasta), 2)

        out_1 = (self.workdir / "out_1").read_text()
        out_2 = (self.workdir / "out_2").read_text()
        self.assertIn(">chr1 description", out_1)
        self.assertIn(">chr2", out_1)
        self.assertNotIn(">scaffold_3", out_1)
        self.assertEqual(out_2, ">scaffold_3\nTTTT\n")

    def test_extract_based_keywords_supports_exact_and_partial_match(self):
        extract_based_keywords(str(self.fasta), "chr2", "True")
        self.assertEqual((self.workdir / "genome.extract.fasta").read_text(), ">chr2\nNNNNAA\n")

        extract_based_keywords(str(self.fasta), "scaffold", None)
        self.assertEqual((self.workdir / "genome.extract.fasta").read_text(), ">scaffold_3\nTTTT\n")

        extract_based_keywords(str(self.fasta), "chr", "False")
        output = (self.workdir / "genome.extract.fasta").read_text()
        self.assertIn(">chr1 description", output)
        self.assertIn(">chr2", output)
        self.assertNotIn(">scaffold_3", output)

    def test_extract_based_length_and_coords(self):
        extract_based_length(str(self.fasta), 6)
        self.assertEqual((self.workdir / "genome.extract.fasta").read_text(), ">chr1 description\nACGTACGTACGT\n")

        extract_based_coords(str(self.fasta), "chr1 description-4-8")
        self.assertEqual((self.workdir / "genome.extract.fasta").read_text(), ">chr1 description:4-8\nTACGT\n")

    def test_func_get_fasta_length_writes_sorted_and_unsorted_lengths(self):
        with contextlib.redirect_stdout(io.StringIO()):
            func_get_fasta_length(str(self.fasta))

        self.assertEqual(
            (self.workdir / "genome.falen").read_text(),
            "chr1 description\t12\nchr2\t6\nscaffold_3\t4\n",
        )
        self.assertEqual(
            (self.workdir / "genome.sort.falen").read_text(),
            "chr1 description\t12\nchr2\t6\nscaffold_3\t4\n",
        )

    def test_genome_change_chr_name_streams_records(self):
        id_map = self.workdir / "ids.tsv"
        id_map.write_text("chr2\tchrB\nscaffold_3\tscaffoldC\n")

        genome_change_chr_name(str(self.fasta), str(id_map))

        self.assertEqual(
            (self.workdir / "genome.ID.change.fasta").read_text(),
            ">chr1 description\nACGTACGTACGT\n>chrB\nNNNNAA\n>scaffoldC\nTTTT\n",
        )

    def test_genome_reverse_some_chr_streams_records(self):
        ids = self.workdir / "reverse.ids"
        ids.write_text("chr2\n")

        genome_reverse_some_chr(str(self.fasta), str(ids))

        self.assertEqual(
            (self.workdir / "genome.rev.genome.fasta").read_text(),
            ">chr1 description\nACGTACGTACGT\n>chr2\nTTNNNN\n>scaffold_3\nTTTT\n",
        )

    def test_genome_karyotype_streams_lengths(self):
        with contextlib.redirect_stdout(io.StringIO()):
            genome_karyotype(str(self.fasta))

        self.assertEqual(
            (self.workdir / "genome.karyotype.txt").read_text(),
            "chr1 description\t1\t12\nchr2\t1\t6\nscaffold_3\t1\t4\n",
        )

    def test_fasta_extract_based_length_writes_full_sequence(self):
        fasta_extract_based_length(str(self.fasta), 6)

        self.assertEqual(
            (self.workdir / "genome.6.fa").read_text(),
            ">chr1 description\nACGTACGTACGT\n",
        )

    def test_func_fasta_compare_streams_new_fasta(self):
        old_fasta = self.workdir / "old.fa"
        new_fasta = self.workdir / "new.fa"
        old_fasta.write_text(">same\nACGT\n>changed\nAAAA\n")
        new_fasta.write_text(">same\nACGT\n>changed\nAAAT\n>missing\nNN\n")

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            func_fasta_compare(str(new_fasta), str(old_fasta))

        self.assertEqual(
            stdout.getvalue(),
            "same is same\n"
            "changed is not same\n"
            "missing not in " + str(old_fasta) + "\n"
            "missing is not same\n",
        )


if __name__ == "__main__":
    unittest.main()
