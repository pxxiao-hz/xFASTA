# xFASTA development notes

## Local setup

Create a virtual environment and install the project in editable mode:

```shell
python3 -m venv .venv
.venv/bin/pip install -e .
```

Run the test suite before and after refactors:

```shell
.venv/bin/python -m unittest discover -s tests
```

Check the command-line entry points:

```shell
./xFASTA --help
.venv/bin/xfasta --help
```

## Refactor boundaries

Keep changes grouped by feature area so work can continue in separate Codex
threads/worktrees:

- `xFASTA.py`: CLI parsing, model-specific argument validation, command dispatch.
- `util.py` and `fasta_fun.py`: FASTA I/O, sequence extraction, split, length statistics.
- `reads_fun.py`: FASTQ I/O, read ID operations, quality summaries, ONT filtering.
- `gap_fun.py`: gap coordinate reporting in GFF and BED.
- `genome_fun.py`, `gff_fun.py`, `telomere_fun.py`: workflows that call external tools.

Prefer small behavior-preserving steps with focused tests. For large FASTA/FASTQ
files, default to streaming iterators instead of loading whole files into
dictionaries. Keep `fasta_read()` only for legacy paths that still need random
access until they are migrated.

## External commands

New external tool integrations should avoid hard-coded local paths. Prefer CLI
options or config values for tool locations, and use `subprocess.run([...],
check=True)` when a shell pipeline is not required. If a shell pipeline is
unavoidable, keep it isolated and covered by a command-construction test.
