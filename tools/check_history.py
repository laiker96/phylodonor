#!/usr/bin/env python3
"""Scan the whole git history for project patterns, not just the checkout.

A clean working tree says nothing about what is reachable in the repository.
Every blob, commit message, branch and tag name ever pushed stays available
to anyone who clones it, so the history is the thing that is published.

Run before the first push, and in CI on every push:

    python tools/check_history.py

Two sets of patterns are applied.

The **sequence rule** is built in and needs no configuration: a run of 25+
characters from the amino-acid alphabet is a protein sequence or a
nucleotide CDS. This is the rule that matters, because a sequence is the
thing that must not be disclosed and a wordlist of names cannot catch one.

**Project patterns** from ``.private-patterns`` or ``PRIVATE_PATTERNS`` are
applied as well when available, and skipped quietly when not. The other
generic rules in ``check_denylist.py`` (absolute paths, credential shapes)
are deliberately not applied here: they match the repository's own tooling,
which that script scans with its own exclusions.

Exit status is 1 if anything matches.

If this ever fails, note that removing the material from the current tree is
not enough: the history has to be rewritten (``git filter-repo``) and, if the
commits were already pushed to a public remote, treat the content as
disclosed regardless.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Reuse the loader so the two scripts cannot disagree about where patterns
# come from or how they are parsed.
sys.path.insert(0, str(REPO / "tools"))
from check_denylist import BUILTIN, load_project_patterns  # noqa: E402

#: Built-in rules meaningful over history. Only the sequence rule: the
#: others match this repository's own source.
_HISTORY_BUILTIN = [(desc, pattern) for desc, pattern in BUILTIN if "sequence" in desc]


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=REPO,
        capture_output=True,
        text=True,
        errors="ignore",
    ).stdout


def main() -> int:
    project = load_project_patterns()
    if not project:
        print(
            "note: no project patterns configured (.private-patterns absent, "
            "PRIVATE_PATTERNS unset); scanning with the built-in sequence rule only",
            file=sys.stderr,
        )

    rules = _HISTORY_BUILTIN + project
    if not rules:  # pragma: no cover - BUILTIN always has the sequence rule
        return 1
    patterns = [(desc, re.compile(p, re.IGNORECASE)) for desc, p in rules]

    # Each source is scanned separately so a hit names where it lives: a blob
    # needs a history rewrite, a ref name only needs the ref deleted.
    sources = {
        "file contents (all commits)": _git("log", "--all", "-p", "--no-color"),
        "commit messages": _git("log", "--all", "--pretty=format:%s%n%b"),
        "file paths": _git("log", "--all", "--pretty=format:", "--name-only"),
        "branch and tag names": _git("for-each-ref", "--format=%(refname)"),
    }

    findings: list[str] = []
    for where, text in sources.items():
        for desc, rx in patterns:
            if rx.search(text):
                findings.append(f"{where}: /{rx.pattern}/ ({desc})")

    if findings:
        print(f"\nBLOCKED: {len(findings)} finding(s) in git history.\n", file=sys.stderr)
        for f in findings:
            print(f"  {f}", file=sys.stderr)
        print(
            "\nThe history is what gets published. Removing the material from "
            "the current tree is not sufficient -- rewrite the history, and if "
            "these commits were already pushed publicly, treat the content as "
            "disclosed.\n",
            file=sys.stderr,
        )
        return 1

    total = sum(len(t) for t in sources.values())
    print(
        f"history clean: {total:,} chars across {len(sources)} sources, "
        f"{len(patterns)} project pattern(s)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
