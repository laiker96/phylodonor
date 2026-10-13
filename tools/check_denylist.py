#!/usr/bin/env python3
"""Refuse to commit a biological sequence or a structure file.

Runs as a pre-commit hook. This is small on purpose.

The pipeline is generic by design -- no sequence, structure or accession
needs to live in the code, and none does. So this is not defending against a
systemic problem; it exists for one accident. A `git push` to a public remote
cannot be undone: a sequence committed and pushed is disclosed permanently,
whatever is done to the history afterwards. That asymmetry is worth a short
check, and not much more.

Deliberately absent: a scan of the git history (the hook prevents the commit,
so there is nothing to find later) and a CI job (CI runs *after* the push, by
which point a public repository has already published the content -- an
alarm, not a prevention).

    pre-commit install                      # do this once
    python tools/check_denylist.py --all    # or check the whole tree by hand
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# A run of 25+ letters drawn only from the amino-acid alphabet is a protein
# sequence. A/C/G/T are a subset, so the same rule catches a nucleotide
# construct. 25 clears the 20-character alphabet constant that
# sequence-handling code legitimately contains.
SEQUENCE = re.compile(r"(?<![A-Za-z])[ACDEFGHIKLMNPQRSTVWY]{25,}(?![A-Za-z])", re.IGNORECASE)

# File types that are always project data, never source.
FORBIDDEN_PATH = re.compile(
    r"\.(?:pdb|cif|mmcif|ent|mol2|sdf|prmtop|inpcrd|rst7|dcd|xtc"
    r"|fa|fasta|faa|fna|seq|sto|a3m|hmm)$",
    re.IGNORECASE,
)

SKIP_SUFFIXES = {".png", ".jpg", ".gif", ".pdf", ".gz", ".zip", ".whl", ".pt"}
MAX_BYTES = 2_000_000


def tracked_files(staged: bool) -> list[Path]:
    cmd = (
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"]
        if staged
        else ["git", "ls-files"]
    )
    try:
        out = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    return [REPO / line for line in out.splitlines() if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="check every tracked file")
    args = parser.parse_args()

    findings: list[str] = []
    self_path = Path(__file__).resolve()

    for path in tracked_files(staged=not args.all):
        rel = path.relative_to(REPO).as_posix() if path.is_relative_to(REPO) else str(path)

        if FORBIDDEN_PATH.search(rel):
            findings.append(f"{rel}: structure or sequence file")
            continue
        if not path.is_file() or path.resolve() == self_path:
            continue
        if path.suffix.lower() in SKIP_SUFFIXES or path.stat().st_size > MAX_BYTES:
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            if SEQUENCE.search(line):
                findings.append(f"{rel}:{i}: looks like a biological sequence")
                break

    if findings:
        print(f"\nBLOCKED: {len(findings)} finding(s).\n", file=sys.stderr)
        for finding in findings:
            print(f"  {finding}", file=sys.stderr)
        print(
            "\nSequences, structures and designed variants belong in your input "
            "directories, which are gitignored -- not in the repository.\n",
            file=sys.stderr,
        )
        return 1

    print("clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())
