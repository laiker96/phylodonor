#!/usr/bin/env python3
"""Pre-publication scan: refuse to commit project-specific or private material.

Run as a pre-commit hook over staged files, and in CI over the whole tree.

    python tools/check_denylist.py --staged
    python tools/check_denylist.py --all

Two classes of pattern are checked.

**Built-in patterns** (below) are generic: personal absolute paths, credential
shapes, REPL history files, redistributed PDFs. They ship with the repo
because they leak nothing.

**Project patterns** are yours: the accessions, organism names, internal
codenames and sequence identifiers that must not reach a public tree. A
pattern list naming the thing it protects is itself a disclosure, so these are
never committed -- keep them in ``.private-patterns``, which is gitignored,
one case-insensitive regex per line, ``#`` for comments:

    # .private-patterns  -- NEVER COMMIT THIS FILE
    <internal-codename>
    <accession>
    <organism-name>

In CI, provide the same content through the ``PRIVATE_PATTERNS`` environment
variable (from a repository secret) rather than a file.

Note that the placeholders above are deliberately not real: project patterns
are applied to this file too, precisely so that nobody documents the scanner
by pasting in the material it is meant to catch.

Exit status is 1 if anything matches, so the commit or build fails.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DENYLIST_FILE = REPO / ".private-patterns"

# Generic patterns, safe to publish. (description, regex)
BUILTIN: list[tuple[str, str]] = [
    ("personal absolute path", r"/(?:Users|home)/(?!runner\b|user\b)[A-Za-z0-9._-]+/"),
    ("Windows user path", r"[A-Z]:\\\\Users\\\\[A-Za-z0-9._-]+"),
    ("private key block", r"-----BEGIN (?:RSA |EC |OPENSSH |PGP )?PRIVATE KEY"),
    ("AWS access key id", r"\bAKIA[0-9A-Z]{16}\b"),
    (
        "generic API token assignment",
        r"(?i)\b(?:api[_-]?key|secret|passwd|password|token)\b\s*[:=]\s*['\"][^'\"\s]{12,}",
    ),
    ("GitHub token", r"\bgh[pousr]_[A-Za-z0-9]{20,}\b"),
    ("bearer token header", r"(?i)authorization\s*[:=]\s*['\"]?bearer\s+\S{12,}"),
    (
        "hardcoded email address",
        r"\b[A-Za-z0-9._%+-]+@(?!example\.(?:com|org)\b)[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b",
    ),
]

# Paths that must never be tracked, regardless of content. Checked before the
# content scan, so a force-added file is still caught.
FORBIDDEN_PATHS: list[tuple[str, str]] = [
    # The pattern file names the material it protects, so committing it is
    # the leak it exists to prevent.
    ("private pattern file (must stay gitignored)", r"(?:^|/)\.private-patterns$"),
    # A tracked file overrides .gitignore, so `git add -f` or `git mv` on one
    # of these silently defeats the ignore rule. Catch it here instead.
    ("private note file (must stay gitignored)", r"(?:^|/)(?:PRIVATE[^/]*|NOTES\.local)\.md$"),
    # Matches both `.Rhistory` and `anything.Rhistory`.
    (
        "shell/REPL history",
        r"(?:^|/|\.)(?:Rhistory|Rapp\.history|bash_history|zsh_history|python_history)$",
    ),
    ("environment file", r"(?:^|/)\.env(?:\.|$)"),
    ("redistributed PDF", r"\.pdf$"),
    ("proprietary instrument data", r"\.(?:dxpx|pse|sdfx)$"),
    ("conda env directory", r"(?:^|/)envs/[^/]+/(?:bin|lib|conda-meta)/"),
]

# Extensions we do not grep for content (binary or large).
SKIP_CONTENT = {
    ".png",
    ".jpg",
    ".jpeg",
    ".gif",
    ".ico",
    ".pdf",
    ".gz",
    ".zip",
    ".bz2",
    ".xz",
    ".so",
    ".dylib",
    ".whl",
    ".pt",
    ".ckpt",
    ".npz",
}

MAX_BYTES = 2_000_000


def load_project_patterns() -> list[tuple[str, str]]:
    """Read project-specific patterns from env var or gitignored file."""
    raw = os.environ.get("PRIVATE_PATTERNS")
    source = "PRIVATE_PATTERNS env var"
    if raw is None and DENYLIST_FILE.exists():
        raw = DENYLIST_FILE.read_text()
        source = str(DENYLIST_FILE.name)
    if not raw:
        return []
    out = []
    for line in raw.splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append((f"project denylist ({source})", line))
    return out


def tracked_files(staged: bool) -> list[Path]:
    cmd = (
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"]
        if staged
        else ["git", "ls-files"]
    )
    try:
        res = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    return [REPO / line for line in res.stdout.splitlines() if line.strip()]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--staged", action="store_true", help="check staged changes only")
    group.add_argument("--all", action="store_true", help="check every tracked file")
    ap.add_argument(
        "--require-project-denylist",
        action="store_true",
        help="fail if no project patterns are configured (use in CI)",
    )
    args = ap.parse_args()

    project = load_project_patterns()
    if args.require_project_denylist and not project:
        print(
            "FAIL: no project denylist configured. Set PRIVATE_PATTERNS or create "
            f"{DENYLIST_FILE.name}.",
            file=sys.stderr,
        )
        return 1
    if not project:
        print(
            f"note: no project denylist found ({DENYLIST_FILE.name} absent, "
            "PRIVATE_PATTERNS unset) -- only generic patterns were checked.",
            file=sys.stderr,
        )

    builtin_rx = [(d, re.compile(p, re.IGNORECASE)) for d, p in BUILTIN]
    project_rx = [(d, re.compile(p, re.IGNORECASE)) for d, p in project]
    path_rules = [(d, re.compile(p, re.IGNORECASE)) for d, p in FORBIDDEN_PATHS]

    files = tracked_files(staged=args.staged or not args.all)
    findings: list[str] = []
    # This file defines the built-in regexes, so they self-match here. Project
    # patterns are still applied to it -- this script has no business
    # containing project identifiers, and documenting it by pasting in real
    # ones is an easy mistake that nothing else would catch.
    self_path = Path(__file__).resolve()

    for path in files:
        rel = path.relative_to(REPO).as_posix() if path.is_relative_to(REPO) else str(path)

        for desc, rx in path_rules:
            if rx.search(rel):
                findings.append(f"{rel}: forbidden path type -- {desc}")

        if not path.is_file():
            continue
        # The pattern file is every pattern by construction; scanning its
        # contents yields nothing but self-matches. Its path rule above is
        # what keeps it out of the tree.
        if path.resolve() == DENYLIST_FILE.resolve():
            continue
        if path.suffix.lower() in SKIP_CONTENT or path.stat().st_size > MAX_BYTES:
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue

        patterns = project_rx if path.resolve() == self_path else builtin_rx + project_rx
        for i, line in enumerate(text.splitlines(), 1):
            for desc, rx in patterns:
                if rx.search(line):
                    findings.append(f"{rel}:{i}: {desc} -- /{rx.pattern}/")
                    break

    if findings:
        print(f"\nBLOCKED: {len(findings)} finding(s).\n", file=sys.stderr)
        for f in findings:
            print(f"  {f}", file=sys.stderr)
        print(
            "\nRemove the material, or narrow the pattern if this is a false positive.\n",
            file=sys.stderr,
        )
        return 1

    print(
        f"clean: {len(files)} file(s) checked, "
        f"{len(builtin_rx)} built-in + {len(project_rx)} project pattern(s), "
        f"{len(path_rules)} path rule(s)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
