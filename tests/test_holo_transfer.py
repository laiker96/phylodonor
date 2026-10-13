"""Tests for cofactor transfer.

The central test is a round trip: a donor built by rigidly moving the
receptor and attaching a ligand must, after transfer, put that ligand back on
its original coordinates. That single assertion covers sequence alignment,
superposition, the coordinate transform and the PDB round trip at once -- if
any of them is wrong the ligand lands somewhere else.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from conftest import write_structure

from phylodonor.holo import fetch, transfer
from phylodonor.holo.structure import parse_atom_records
from phylodonor.holo.transfer import AnnotationSpec, CofactorSpec


def _spec(donor: Path, **kwargs) -> CofactorSpec:
    return CofactorSpec(
        resname="LIG",
        donor_pdb="1ABC",
        explicit_path=donor,
        min_fit_pairs=40,
        **kwargs,
    )


def test_transfer_recovers_ligand_position(tmp_path, receptor_pdb, donor_pdb, ligand_coords):
    """A rigid-transform donor must place its ligand back at the original spot."""
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb)],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )

    atoms = [a for a in parse_atom_records(result.holo_pdb) if a.resn == "LIG"]
    assert len(atoms) == len(ligand_coords)

    placed = np.array([a.coord for a in atoms])
    # PDB coordinates carry three decimals, so ~1e-3 is the floor here.
    assert np.allclose(placed, ligand_coords, atol=5e-3), (
        f"max deviation {np.abs(placed - ligand_coords).max():.4f} A"
    )


def test_transfer_fit_is_near_perfect(tmp_path, receptor_pdb, donor_pdb):
    """Identical geometry under a rigid transform should fit to ~0 RMSD."""
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb)],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )
    assert result.fits["lig_global_rmsd_a"] < 0.01
    assert result.fits["lig_global_pairs"] == 60


def test_protein_is_written_on_configured_chain(tmp_path, receptor_pdb, donor_pdb):
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb)],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )
    atoms = parse_atom_records(result.holo_pdb)
    assert {a.chain for a in atoms if a.kind == "ATOM"} == {"A"}
    # The cofactor must not collide with the protein chain.
    assert result.cofactor_chains["LIG"] == "B"
    assert {a.chain for a in atoms if a.resn == "LIG"} == {"B"}


def test_serials_are_contiguous_and_ter_separates_ligand(tmp_path, receptor_pdb, donor_pdb):
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb)],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )
    lines = result.holo_pdb.read_text().splitlines()
    serials = [int(line[6:11]) for line in lines if line.startswith(("ATOM", "HETATM"))]
    assert serials == list(range(1, len(serials) + 1))
    assert "TER" in lines
    assert lines.index("TER") < next(i for i, line in enumerate(lines) if line.startswith("HETATM"))


def test_explicit_out_chain_is_honoured(tmp_path, receptor_pdb, donor_pdb):
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb, out_chain="Z")],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )
    assert result.cofactor_chains["LIG"] == "Z"


def test_two_cofactors_get_distinct_chains(tmp_path, receptor_pdb, ligand_coords, donor_transform):
    donor = write_structure(
        tmp_path / "donor_two.pdb",
        n_residues=60,
        protein_chain="X",
        ligands={
            "Y": [("LIG", c, "C") for c in ligand_coords],
            "Z": [("CO2", c + np.array([0.0, 8.0, 0.0]), "C") for c in ligand_coords],
        },
        transform=donor_transform,
    )
    specs = [
        CofactorSpec(resname="LIG", donor_pdb="1ABC", explicit_path=donor, min_fit_pairs=40),
        CofactorSpec(resname="CO2", donor_pdb="1ABC", explicit_path=donor, min_fit_pairs=40),
    ]
    result = transfer.transfer_cofactors(
        receptor_pdb, specs, [], tmp_path / "out", model_id="demo", allow_download=False
    )
    assert result.cofactor_chains == {"LIG": "B", "CO2": "C"}
    assert len({*result.cofactor_chains.values()}) == 2


def test_annotation_ligands_go_only_in_the_annotation_model(
    tmp_path, receptor_pdb, ligand_coords, donor_transform
):
    donor = write_structure(
        tmp_path / "donor_anno.pdb",
        n_residues=60,
        protein_chain="X",
        ligands={
            "Y": [("LIG", c, "C") for c in ligand_coords],
            "W": [("ANN", c + np.array([0.0, 0.0, 9.0]), "N") for c in ligand_coords],
        },
        transform=donor_transform,
    )
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor)],
        [AnnotationSpec(label="substrate", resnames=("ANN",), from_donor="LIG")],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )

    holo_resns = {a.resn for a in parse_atom_records(result.holo_pdb)}
    assert "LIG" in holo_resns
    assert "ANN" not in holo_resns, "annotation ligands must never reach the minimised model"

    assert result.annotation_pdb is not None
    anno_resns = {a.resn for a in parse_atom_records(result.annotation_pdb)}
    assert {"LIG", "ANN"} <= anno_resns
    assert result.annotation_counts["substrate"] == len(ligand_coords)


def test_no_annotation_model_when_none_configured(tmp_path, receptor_pdb, donor_pdb):
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb)],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )
    assert result.annotation_pdb is None


def test_missing_cofactor_in_donor_is_an_error(tmp_path, receptor_pdb, donor_pdb):
    spec = CofactorSpec(resname="XYZ", donor_pdb="1ABC", explicit_path=donor_pdb, min_fit_pairs=40)
    with pytest.raises(ValueError, match="does not contain ligand"):
        transfer.transfer_cofactors(
            receptor_pdb,
            [spec],
            [],
            tmp_path / "out",
            model_id="demo",
            allow_download=False,
        )


def test_too_few_pairs_refuses_rather_than_guessing(tmp_path, ligand_coords, donor_transform):
    """A fit on a handful of residues is arbitrary and must not silently pass."""
    short_receptor = write_structure(tmp_path / "short.pdb", n_residues=20)
    donor = write_structure(
        tmp_path / "donor_short.pdb",
        n_residues=20,
        protein_chain="X",
        ligands={"Y": [("LIG", c, "C") for c in ligand_coords]},
        transform=donor_transform,
    )
    with pytest.raises(ValueError, match="too few aligned CA pairs"):
        transfer.transfer_cofactors(
            short_receptor,
            [_spec(donor)],
            [],
            tmp_path / "out",
            model_id="demo",
            allow_download=False,
        )


def test_missing_input_chain_is_an_error(tmp_path, receptor_pdb, donor_pdb):
    with pytest.raises(ValueError, match="not found in"):
        transfer.transfer_cofactors(
            receptor_pdb,
            [_spec(donor_pdb)],
            [],
            tmp_path / "out",
            model_id="demo",
            input_chain="Q",
            allow_download=False,
        )


def test_empty_cofactor_list_is_an_error(tmp_path, receptor_pdb):
    with pytest.raises(ValueError, match="no cofactors configured"):
        transfer.transfer_cofactors(
            receptor_pdb, [], [], tmp_path / "out", model_id="demo", allow_download=False
        )


def test_downloads_disabled_reports_clearly(tmp_path, receptor_pdb):
    spec = CofactorSpec(resname="LIG", donor_pdb="9ZZZ")
    with pytest.raises(fetch.DonorUnavailable, match="downloads are disabled"):
        transfer.transfer_cofactors(
            receptor_pdb,
            [spec],
            [],
            tmp_path / "out",
            model_id="demo",
            allow_download=False,
        )


def test_cached_donor_is_used_without_network(tmp_path, receptor_pdb, donor_pdb):
    """A donor pre-staged in the cache must be found by pdb id alone."""
    out_dir = tmp_path / "out"
    cache = out_dir / "donor_crystals"
    cache.mkdir(parents=True)
    (cache / "1ABC.pdb").write_text(donor_pdb.read_text())

    result = transfer.transfer_cofactors(
        receptor_pdb,
        [CofactorSpec(resname="LIG", donor_pdb="1ABC", min_fit_pairs=40)],
        [],
        out_dir,
        model_id="demo",
        allow_download=False,
    )
    assert result.donors["LIG"]["source"] == "cached"


def test_pocket_refit_runs_and_is_reported(tmp_path, receptor_pdb, donor_pdb):
    result = transfer.transfer_cofactors(
        receptor_pdb,
        [_spec(donor_pdb, pocket_cutoff_angstrom=8.0)],
        [],
        tmp_path / "out",
        model_id="demo",
        allow_download=False,
    )
    assert "lig_pocket_rmsd_a" in result.fits
    assert result.fits["lig_pocket_pairs"] >= 4
    # The pocket is a subset of the protein, so it fits on fewer pairs.
    assert result.fits["lig_pocket_pairs_considered"] < result.fits["lig_global_pairs_considered"]
