import contextlib
import io
import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import xFASTA  # noqa: E402


class CliTest(unittest.TestCase):
    def parse_with_stderr(self, argv):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as error:
                xFASTA.main(argv)
        return error.exception.code, stderr.getvalue()

    def test_importing_cli_does_not_load_heavy_bio_modules(self):
        code = (
            "import sys; import xFASTA; "
            "print('Bio' in sys.modules, 'pandas' in sys.modules, 'pysam' in sys.modules)"
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            capture_output=True,
            check=True,
        )

        self.assertEqual(result.stdout.strip(), "False False False")

    def test_help_returns_without_dispatching(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            with self.assertRaises(SystemExit) as error:
                xFASTA.main(["--help"])

        self.assertEqual(error.exception.code, 0)
        self.assertIn("处理 FASTA 文件", stdout.getvalue())

    def test_fasta_model_is_required(self):
        code, stderr = self.parse_with_stderr(["fasta", "-i", "genome.fa"])

        self.assertEqual(code, 2)
        self.assertIn("the following arguments are required: -m/--model", stderr)

    def test_fasta_extract_requires_keywords(self):
        code, stderr = self.parse_with_stderr(["fasta", "-i", "genome.fa", "-m", "extractk"])

        self.assertEqual(code, 2)
        self.assertIn("fasta extractk requires -k/--keywords", stderr)

    def test_fasta_extractk_defaults_to_partial_keyword_match(self):
        with mock.patch("fasta_fun.extract_based_keywords") as extract:
            exit_code = xFASTA.main(["fasta", "-i", "genome.fa", "-m", "extractk", "-k", "chr"])

        self.assertEqual(exit_code, 0)
        extract.assert_called_once_with("genome.fa", "chr", "False")

    def test_read_extract_requires_read_id(self):
        code, stderr = self.parse_with_stderr(["read", "-i", "reads.fastq.gz", "-m", "extract"])

        self.assertEqual(code, 2)
        self.assertIn("read extract requires -k/--readid", stderr)

    def test_read_filteront_dispatches_configured_chopper(self):
        with mock.patch("reads_fun.func_reads_filter_ONT_reads") as filter_reads:
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = xFASTA.main([
                    "read",
                    "-i",
                    "reads.fastq.gz",
                    "-m",
                    "filterONT",
                    "-q",
                    "7",
                    "-l",
                    "100000",
                    "--chopper",
                    "/opt/chopper",
                ])

        self.assertEqual(exit_code, 0)
        filter_reads.assert_called_once_with("reads.fastq.gz", 7, 100000, "/opt/chopper")

    def test_genome_split_by_gap_requires_input_genome(self):
        code, stderr = self.parse_with_stderr(["genome", "-m", "split_by_gap"])

        self.assertEqual(code, 2)
        self.assertIn("genome split_by_gap requires -i/--inputGenome", stderr)

    def test_telomere_dispatches_configured_script(self):
        with mock.patch("telomere_fun.func_telomere_info") as telomere_info:
            exit_code = xFASTA.main([
                "telomere",
                "-i",
                "genome.fa",
                "-o",
                "telo.tsv",
                "--telomere-script",
                "/opt/find_telomere.py",
            ])

        self.assertEqual(exit_code, 0)
        telomere_info.assert_called_once_with(
            "genome.fa", 150000, 100, "telo.tsv", "/opt/find_telomere.py"
        )

    def test_genome_telomere_dispatches_configured_script(self):
        with mock.patch("genome_fun.genome_telomere") as genome_telomere:
            exit_code = xFASTA.main([
                "genome",
                "-m",
                "telomere",
                "-i",
                "genome.fa",
                "--telomere-script",
                "/opt/find_telomere.py",
            ])

        self.assertEqual(exit_code, 0)
        genome_telomere.assert_called_once_with("genome.fa", "/opt/find_telomere.py")

    def test_genome_survey_dispatches_configured_genomescope_scripts(self):
        with mock.patch("genome_fun.genome_survey") as genome_survey:
            exit_code = xFASTA.main([
                "genome",
                "-m",
                "survey",
                "-p",
                "reads",
                "--genomescope1-script",
                "/opt/gs1.R",
                "--genomescope2-script",
                "/opt/gs2.R",
            ])

        self.assertEqual(exit_code, 0)
        genome_survey.assert_called_once_with(
            21,
            150,
            10000000,
            16,
            "reads",
            2,
            "/opt/gs1.R",
            "/opt/gs2.R",
        )


if __name__ == "__main__":
    unittest.main()
