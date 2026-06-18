import subprocess
import sys
import unittest
from unittest import mock

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from telomere_fun import build_telomere_command, func_telomere_info  # noqa: E402


class TelomereFunTest(unittest.TestCase):
    def test_build_telomere_command_uses_configured_script(self):
        self.assertEqual(
            build_telomere_command("genome.fa", 150000, 100, "telo.tsv", "/opt/find_telomere.py", "python3"),
            [
                "python3",
                "/opt/find_telomere.py",
                "-i", "genome.fa",
                "-l", "150000",
                "-m", "100",
                "-o", "telo.tsv",
            ],
        )

    def test_func_telomere_info_uses_subprocess_without_shell(self):
        with mock.patch("telomere_fun.subprocess.run") as run:
            func_telomere_info("genome.fa", 150000, 100, "telo.tsv", "/opt/find_telomere.py")

        command = run.call_args.args[0]
        self.assertEqual(command[1:], [
            "/opt/find_telomere.py",
            "-i", "genome.fa",
            "-l", "150000",
            "-m", "100",
            "-o", "telo.tsv",
        ])
        self.assertEqual(run.call_args.kwargs, {"close_fds": True, "check": True})


if __name__ == "__main__":
    unittest.main()
