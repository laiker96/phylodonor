"""Tests for contact annotation.

The interesting behaviour is the clash heuristic, which has to distinguish a
hydrogen bond from a steric problem. Both show a large formal van der Waals
overlap; only one is a defect. Getting this backwards would condemn every
correctly placed cofactor, so it is tested directly rather than through the
pipeline.
"""

from __future__ import annotations

import numpy as np
import pytest

from phylodonor.holo.contacts import (
    FIELDS,
    LigandGroup,
    closest_approach,
    contact_rows,
    write_contacts,
)
from phylodonor.holo.structure import AtomRecord


def atom(name, elem, coord, *, kind="ATOM", resn="ALA", chain="A", resi="1"):
    return AtomRecord(
        kind=kind,
        name=name,
        resn=resn,
        chain=chain,
        resi=resi,
        icode="",
        elem=elem,
        coord=np.array(coord, dtype=float),
    )


def test_closest_approach_reports_nearest_pair():
    protein = [atom("CA", "C", [0, 0, 0]), atom("CB", "C", [10, 0, 0])]
    target = [atom("C1", "C", [3.5, 0, 0], kind="HETATM", resn="LIG")]
    geom = closest_approach(protein, target)
    assert geom is not None
    assert geom["min_distance_a"] == pytest.approx(3.5)
    assert geom["protein_atom"] == "CA"


def test_hydrogen_bond_is_not_a_severe_clash():
    """N and O at 2.7 A overlap formally but are a textbook hydrogen bond."""
    protein = [atom("ND1", "N", [0, 0, 0], resn="HIS")]
    target = [atom("O2", "O", [2.7, 0, 0], kind="HETATM", resn="LIG")]
    geom = closest_approach(protein, target)
    assert geom is not None
    # 1.55 + 1.52 - 2.7 = 0.37, below the 0.4 reporting threshold.
    assert geom["severe_clash"] is False


def test_close_polar_contact_is_flagged_hbond_compatible():
    protein = [atom("ND1", "N", [0, 0, 0], resn="HIS")]
    target = [atom("O2", "O", [2.4, 0, 0], kind="HETATM", resn="LIG")]
    geom = closest_approach(protein, target)
    assert geom is not None
    assert geom["max_vdw_overlap_a"] == pytest.approx(1.55 + 1.52 - 2.4)
    assert geom["hbond_compatible"] is True
    assert geom["severe_clash"] is False


def test_carbon_overlap_is_a_severe_clash():
    """Two carbons far inside each other's radii cannot be explained away."""
    protein = [atom("CB", "C", [0, 0, 0])]
    target = [atom("C1", "C", [2.6, 0, 0], kind="HETATM", resn="LIG")]
    geom = closest_approach(protein, target)
    assert geom is not None
    assert geom["max_vdw_overlap_a"] == pytest.approx(1.70 + 1.70 - 2.6)
    assert geom["severe_clash"] is True
    assert geom["hbond_compatible"] is False


def test_polar_pair_too_close_to_be_an_hbond_is_severe():
    """Below 2.2 A even N/O is not a hydrogen bond, it is a collision."""
    protein = [atom("ND1", "N", [0, 0, 0], resn="HIS")]
    target = [atom("O2", "O", [1.5, 0, 0], kind="HETATM", resn="LIG")]
    geom = closest_approach(protein, target)
    assert geom is not None
    assert geom["severe_clash"] is True
    assert geom["hbond_compatible"] is False


def test_hydrogens_are_ignored():
    protein = [atom("HA", "H", [1.0, 0, 0]), atom("CA", "C", [5.0, 0, 0])]
    target = [atom("C1", "C", [0, 0, 0], kind="HETATM", resn="LIG")]
    geom = closest_approach(protein, target)
    assert geom is not None
    assert geom["protein_atom"] == "CA", "hydrogen should not win the closest pair"


def test_empty_side_returns_none():
    assert closest_approach([], [atom("C1", "C", [0, 0, 0])]) is None
    assert closest_approach([atom("CA", "C", [0, 0, 0])], []) is None


def test_contact_rows_respects_cutoff_and_sorts_by_residue():
    protein = [
        atom("CA", "C", [0, 0, 0], resi="7"),
        atom("CA", "C", [3.0, 0, 0], resi="3"),
        atom("CA", "C", [50.0, 0, 0], resi="5"),
    ]
    group = LigandGroup(
        label="lig_contacts",
        atoms=(atom("C1", "C", [1.0, 0, 0], kind="HETATM", resn="LIG", chain="B", resi="1"),),
    )
    rows = contact_rows(protein, group, cutoff_a=4.0)
    assert [r.resi for r in rows] == ["3", "7"], "far residue excluded, rest sorted numerically"
    assert all(r.annotation == "lig_contacts" for r in rows)
    assert rows[0].target_chain == "B"


def test_contact_rows_only_considers_the_protein_chain():
    protein = [
        atom("CA", "C", [0, 0, 0], chain="A", resi="1"),
        atom("CA", "C", [0.5, 0, 0], chain="B", resi="2"),
    ]
    group = LigandGroup(label="x", atoms=(atom("C1", "C", [1.0, 0, 0], kind="HETATM", resn="LIG"),))
    rows = contact_rows(protein, group, cutoff_a=4.0, protein_chain="A")
    assert [r.chain for r in rows] == ["A"]


def test_hetatm_protein_records_are_not_counted_as_protein():
    protein = [atom("C1", "C", [0, 0, 0], kind="HETATM", resn="LIG")]
    group = LigandGroup(label="x", atoms=(atom("C2", "C", [1.0, 0, 0], kind="HETATM", resn="LIG"),))
    assert contact_rows(protein, group, cutoff_a=4.0) == []


def test_write_contacts_emits_header_even_when_empty(tmp_path):
    out = write_contacts([], tmp_path / "contacts.tsv")
    lines = out.read_text().splitlines()
    assert lines == ["\t".join(FIELDS)], "an empty result is still a result"


def test_write_contacts_round_trips(tmp_path):
    protein = [atom("CA", "C", [0, 0, 0], resi="12", resn="TRP")]
    group = LigandGroup(
        label="lig_contacts",
        atoms=(atom("C1", "C", [3.0, 0, 0], kind="HETATM", resn="LIG", chain="B", resi="1"),),
    )
    rows = contact_rows(protein, group, cutoff_a=4.0)
    out = write_contacts(rows, tmp_path / "contacts.tsv")

    header, body = out.read_text().splitlines()
    assert header.split("\t") == list(FIELDS)
    record = dict(zip(FIELDS, body.split("\t"), strict=True))
    assert record["resi"] == "12"
    assert record["resn"] == "TRP"
    assert record["min_distance_a"] == "3.000"
    assert record["severe_clash"] == "no"
