import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from util import check_file_in_path, func_n50, get_Nnumber, rev_seq  # noqa: E402


class UtilTest(unittest.TestCase):
    def test_check_file_in_path_skips_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            existing = Path(tmp) / "done.txt"
            existing.write_text("ok")

            with mock.patch("util.subprocess.run") as run:
                with contextlib.redirect_stdout(io.StringIO()):
                    check_file_in_path(str(existing), "touch done.txt")

        run.assert_not_called()

    def test_check_file_in_path_checks_command_return_code(self):
        with mock.patch("util.subprocess.run") as run:
            check_file_in_path("missing.txt", "false")

        run.assert_called_once_with("false", shell=True, check=True)

    def test_get_nnumber_does_not_mutate_lengths(self):
        lengths = [4, 10, 6]

        self.assertEqual(get_Nnumber(20, lengths, 50), (10, 1))
        self.assertEqual(lengths, [4, 10, 6])

    def test_func_n50_rejects_empty_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            previous_cwd = Path.cwd()
            try:
                os.chdir(tmp)
                empty = Path(tmp) / "empty.fa"
                empty.write_text("")

                with self.assertRaisesRegex(ValueError, "no sequences"):
                    func_n50(str(empty))
            finally:
                os.chdir(previous_cwd)

    def test_rev_seq_supports_iupac_bases(self):
        self.assertEqual(rev_seq("ACGTRYMKBDHVN"), "NBDHVMKRYACGT")


if __name__ == "__main__":
    unittest.main()
