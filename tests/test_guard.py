"""Tests for the pre-publication guard.

The one rule that matters is the sequence detector: the thing that must not
be published is the target's sequence, its structure, or a designed variant,
and a list of names cannot catch a pasted sequence.

No sequence literal appears in this file, real or synthetic. The positive
cases are *built* from the alphabet at import time, so the guard can scan
its own tests like any other file -- exempting them would mean a real
sequence pasted into a test could slip past, which is the wrong trade for
the sake of a readable fixture.
"""

from __future__ import annotations

import random
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from check_denylist import BUILTIN, FORBIDDEN_PATHS  # noqa: E402


def rule(fragment: str) -> re.Pattern[str]:
    """The compiled built-in rule whose description contains ``fragment``."""
    matches = [p for desc, p in BUILTIN if fragment in desc]
    assert len(matches) == 1, f"expected exactly one rule matching {fragment!r}"
    return re.compile(matches[0], re.IGNORECASE)


def path_rule(fragment: str) -> re.Pattern[str]:
    matches = [p for desc, p in FORBIDDEN_PATHS if fragment in desc]
    assert len(matches) == 1, f"expected exactly one path rule matching {fragment!r}"
    return re.compile(matches[0], re.IGNORECASE)


# --- the sequence detector -------------------------------------------------

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
NUCLEOTIDES = "ACGT"


def synthetic(alphabet: str, length: int, *, seed: int = 0) -> str:
    """A reproducible string over ``alphabet``, built rather than written."""
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(length))


SEQUENCE_LEAKS = [
    # The shape of a pasted FASTA record body.
    synthetic(ALPHABET, 55, seed=1),
    # A nucleotide CDS: A/C/G/T are a subset of the same alphabet, so the one
    # rule covers codon-optimised constructs too.
    synthetic(NUCLEOTIDES, 54, seed=2),
    # Exactly at the 25-character threshold.
    ALPHABET + ALPHABET[:5],
    # Lowercase, as soft-masked sequence files are written.
    synthetic(ALPHABET, 32, seed=3).lower(),
]


@pytest.mark.parametrize("sequence", SEQUENCE_LEAKS)
def test_sequence_detector_catches_sequences(sequence):
    assert rule("biological sequence").search(sequence), sequence


LEGITIMATE = [
    # The 20-character alphabet constant that sequence-handling code uses.
    # Deliberately below the threshold, because this must not be flagged.
    f'AA_ALPHABET = "{ALPHABET}"',
    'SHA = "a3f9c1e8b7d2049fbe15c3a7d9082b4615cfae73"',
    "max_backbone_rmsd_angstrom: 1.5",
    "forbid_n_glycosylation_motifs: false",
    'label = f"{resname.lower()}_contacts"',
    "STANDARD_PROTEIN_RESNAMES = frozenset({'ALA', 'ARG', 'ASN'})",
    "restraint_k_kcal_per_a2: 100.0",
]


@pytest.mark.parametrize("line", LEGITIMATE)
def test_sequence_detector_leaves_ordinary_code_alone(line):
    assert not rule("biological sequence").search(line), line


def test_threshold_is_above_the_alphabet_constant():
    """The 20-letter alphabet is legitimate; 25 is the floor deliberately."""
    pattern = rule("biological sequence")
    assert not pattern.search(ALPHABET)
    assert pattern.search(ALPHABET + ALPHABET[:5])


def test_sequence_embedded_in_a_longer_token_is_still_caught():
    assert rule("biological sequence").search('seq = "' + "A" * 40 + '"')


# --- path rules ------------------------------------------------------------


@pytest.mark.parametrize(
    "path",
    [
        "model.pdb",
        "fold.cif",
        "x.mmcif",
        "ligand.mol2",
        "system.prmtop",
        "traj.dcd",
        "results/run/final.pdb",
    ],
)
def test_structure_files_are_refused(path):
    assert path_rule("structure coordinates").search(path), path


@pytest.mark.parametrize(
    "path",
    ["target.fa", "hits.fasta", "x.faa", "aln.a3m", "profile.hmm", "tree.nwk", "t.treefile"],
)
def test_sequence_files_are_refused(path):
    assert path_rule("biological sequence file").search(path), path


@pytest.mark.parametrize(
    "path",
    ["src/phylodonor/holo/structure.py", "docs/data_sources.md", "configs/template.yaml"],
)
def test_source_files_are_not_refused(path):
    for desc in ("structure coordinates", "biological sequence file"):
        assert not path_rule(desc).search(path), f"{path} matched {desc}"


# --- the repository's own cleanliness --------------------------------------


def test_this_repository_has_no_sequences_or_structures():
    """A regression test on the published tree itself.

    If this fails, something committed a sequence or a structure file, and it
    needs removing from history rather than just from the tree.
    """
    import subprocess

    tracked = subprocess.run(
        ["git", "ls-files"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()

    sequence_rule = rule("biological sequence")
    structures = path_rule("structure coordinates")
    sequences = path_rule("biological sequence file")

    offenders = []
    for name in tracked:
        if structures.search(name) or sequences.search(name):
            offenders.append(f"{name}: forbidden file type")
            continue
        path = REPO_ROOT / name
        if not path.is_file() or path.suffix in {".png", ".gz", ".whl"}:
            continue
        for i, line in enumerate(path.read_text(errors="ignore").splitlines(), 1):
            if sequence_rule.search(line):
                offenders.append(f"{name}:{i}")
    assert not offenders, "sequence or structure material in the tree: " + ", ".join(offenders)
