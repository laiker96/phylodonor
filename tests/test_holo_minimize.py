"""Tests for restrained minimisation.

Skipped unless the ``dynamics`` extra is installed. CI installs it, so these
do run there -- an optional dependency is not a reason to leave the code
untested.

Note on the numbers: the synthetic helix has idealised, slightly unphysical
bond geometry, so the absolute energies and displacements are much larger
than a real predicted structure would give. The *relationships* between them
are what these tests assert, and those hold regardless.
"""

from __future__ import annotations

import numpy as np
import pytest
from conftest import helix_ca_coords, write_structure

from phylodonor.holo import transfer
from phylodonor.holo.structure import parse_atom_records
from phylodonor.holo.transfer import CofactorSpec

pytest.importorskip("openmm", reason="needs the 'dynamics' extra")

from phylodonor.holo import minimize  # noqa: E402

N_RESIDUES = 40


@pytest.fixture
def holo_model(tmp_path):
    """A transferred holo model with a cofactor placed deliberately close."""
    receptor = write_structure(tmp_path / "rec.pdb", n_residues=N_RESIDUES, sequence=["GLY"])
    centre = helix_ca_coords(N_RESIDUES)[N_RESIDUES // 2]
    coords = np.array([centre + [3.0, 0, 0], centre + [4.4, 0, 0]])
    donor = write_structure(
        tmp_path / "don.pdb",
        n_residues=N_RESIDUES,
        sequence=["GLY"],
        protein_chain="X",
        ligands={"Y": [("LIG", c, "C") for c in coords]},
    )
    result = transfer.transfer_cofactors(
        receptor,
        [CofactorSpec(resname="LIG", donor_pdb="1ABC", explicit_path=donor, min_fit_pairs=30)],
        [],
        tmp_path / "out",
        model_id="m",
        allow_download=False,
    )
    return result.holo_pdb


def test_minimisation_lowers_the_energy(holo_model, tmp_path):
    result = minimize.minimize_holo(holo_model, tmp_path / "min.pdb", max_iterations=300)
    assert result.final_energy_kj < result.initial_energy_kj
    assert result.minimized_pdb.exists()


def test_cofactor_is_held_exactly_fixed(holo_model, tmp_path):
    """The cofactor's crystal geometry must survive untouched."""
    before = {(a.resi, a.name): a.coord for a in parse_atom_records(holo_model) if a.resn == "LIG"}
    result = minimize.minimize_holo(holo_model, tmp_path / "min.pdb", max_iterations=300)
    after = {
        (a.resi, a.name): a.coord
        for a in parse_atom_records(result.minimized_pdb)
        if a.resn == "LIG"
    }
    assert before and after.keys() == before.keys()
    for key, coord in before.items():
        assert np.allclose(after[key], coord), f"cofactor atom {key} moved"
    assert result.n_fixed_ligand_atoms == len(before)


def test_stiffer_restraint_moves_the_backbone_less(holo_model, tmp_path):
    """Guards the kcal/mol/A^2 -> kJ/mol/nm^2 conversion.

    A factor-of-418 error there is invisible except as a backbone that
    drifts, so the monotonic relationship is worth pinning.
    """
    soft = minimize.minimize_holo(
        holo_model, tmp_path / "soft.pdb", max_iterations=300, restraint_k_kcal_per_a2=100.0
    )
    stiff = minimize.minimize_holo(
        holo_model, tmp_path / "stiff.pdb", max_iterations=300, restraint_k_kcal_per_a2=1000.0
    )
    free = minimize.minimize_holo(
        holo_model, tmp_path / "free.pdb", max_iterations=300, restrain_backbone=False
    )
    assert stiff.max_ca_displacement_a < soft.max_ca_displacement_a < free.max_ca_displacement_a


def test_atom_count_and_chains_are_preserved(holo_model, tmp_path):
    """Minimisation must not add, drop or relabel atoms."""
    before = parse_atom_records(holo_model)
    result = minimize.minimize_holo(holo_model, tmp_path / "min.pdb", max_iterations=100)
    after = parse_atom_records(result.minimized_pdb)

    # The terminal OXT repair adds at most one atom to the protein.
    assert len(after) - len(before) in (0, 1)
    assert {a.chain for a in after} == {a.chain for a in before}
    assert [a.name for a in after if a.resn == "LIG"] == [a.name for a in before if a.resn == "LIG"]
    assert not any(a.is_hydrogen for a in after), "added hydrogens should not be written out"


def test_serials_are_renumbered_contiguously(holo_model, tmp_path):
    result = minimize.minimize_holo(holo_model, tmp_path / "min.pdb", max_iterations=100)
    lines = result.minimized_pdb.read_text().splitlines()
    serials = [int(line[6:11]) for line in lines if line.startswith(("ATOM", "HETATM"))]
    assert serials == list(range(1, len(serials) + 1))


def test_terminal_oxt_is_added_when_missing(holo_model, tmp_path):
    """Without this repair the force field cannot match the terminal template."""
    assert not any(a.name == "OXT" for a in parse_atom_records(holo_model))
    result = minimize.minimize_holo(holo_model, tmp_path / "min.pdb", max_iterations=50)
    assert any(a.name == "OXT" for a in parse_atom_records(result.minimized_pdb))


def test_no_protein_chain_is_reported_clearly(tmp_path):
    lonely = tmp_path / "ligand_only.pdb"
    lonely.write_text(
        "HETATM    1  C1  LIG B   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n"
    )
    with pytest.raises(minimize.MinimizationUnavailable, match="no ATOM records on chain"):
        minimize.minimize_holo(lonely, tmp_path / "min.pdb")


def test_nonstandard_residue_is_reported_clearly(tmp_path):
    """The error should name the cause, not leak a force-field traceback."""
    weird = write_structure(tmp_path / "weird.pdb", n_residues=30, sequence=["XYZ"])
    with pytest.raises(minimize.MinimizationUnavailable, match="could not parameterise"):
        minimize.minimize_holo(weird, tmp_path / "min.pdb")


def test_manifest_records_the_caveat(holo_model, tmp_path):
    """The output must say what the numbers do and do not mean."""
    result = minimize.minimize_holo(holo_model, tmp_path / "min.pdb", max_iterations=50)
    manifest = result.as_manifest()
    assert "fixed" in manifest["note"] and "uncharged" in manifest["note"]
    assert set(manifest) >= {
        "minimized_pdb",
        "initial_energy_kj_per_mol",
        "final_energy_kj_per_mol",
        "max_ca_displacement_a",
    }
