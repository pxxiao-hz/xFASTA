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

from util import check_file_in_path  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
