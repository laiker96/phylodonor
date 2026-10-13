"""End-to-end tests for the holo stage, driven by a config.

The important one is ``test_pipeline_without_minimisation``: the stage must
produce a usable contact table with minimisation disabled and OpenMM absent,
because that is the configuration the pipeline is meant to run in by default.
"""

from __future__ import annotations

import json
import textwrap

import numpy as np
import pytest
from conftest import helix_ca_coords, write_structure

from phylodonor import config as config_mod
from phylodonor.holo import run as holo_run

N_RESIDUES = 40


def make_config(tmp_path, *, minimize: bool, cutoff: float = 4.5, extra: str = "") -> object:
    """A config pointing at a synthetic receptor and a pre-staged donor."""
    receptor = write_structure(tmp_path / "receptor.pdb", n_residues=N_RESIDUES, sequence=["GLY"])

    centre = helix_ca_coords(N_RESIDUES)[N_RESIDUES // 2]
    coords = np.array([centre + [3.6, 0, 0], centre + [5.0, 0, 0]])
    donor = write_structure(
        tmp_path / "donor.pdb",
        n_residues=N_RESIDUES,
        sequence=["GLY"],
        protein_chain="X",
        ligands={"Y": [("LIG", c, "C") for c in coords]},
    )

    # Pre-stage the donor in the cache the stage will look in, keyed by id,
    # so nothing touches the network.
    cache = tmp_path / "results" / "demo" / "holo" / "receptor" / "donor_crystals"
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "1ABC.pdb").write_text(donor.read_text())

    path = tmp_path / "cfg.yaml"
    path.write_text(
        textwrap.dedent(f"""
        project: {{name: demo}}
        target:
          accession: P00000
          structure: {receptor}
        holo:
          protein_chain: A
          cofactors:
            - resname: LIG
              donor_pdb: 1ABC
              min_fit_pairs: 30
          contacts:
            cutoff_angstrom: {cutoff}
          minimize:
            enabled: {str(minimize).lower()}
            max_iterations: 100
        paths:
          results_root: {tmp_path / "results"}
        {extra}
        """)
    )
    return config_mod.load(path, repo_root=tmp_path)


def test_pipeline_without_minimisation(tmp_path):
    """The default path: transfer and annotate, no OpenMM involved."""
    cfg = make_config(tmp_path, minimize=False)

    built = holo_run.build(cfg, allow_download=False)
    assert built["model_id"] == "receptor"
    out_dir = built["out_dir"]
    assert (out_dir / "receptor_holo.pdb").exists()

    skipped = holo_run.minimize_model(cfg)
    assert skipped["minimized"] is None
    assert not (out_dir / "receptor_minimized.pdb").exists()

    annotated = holo_run.annotate_contacts(cfg)
    assert annotated["contacts_tsv"].exists()
    rows = annotated["rows"]
    assert rows, "the cofactor sits beside the helix, so something must be in contact"
    assert {r.annotation for r in rows} == {"lig_contacts"}

    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert manifest["minimize"] == {"enabled": False}
    assert manifest["contacts"]["minimized"] is False
    assert manifest["contacts"]["model_used"].endswith("receptor_holo.pdb")
    assert "phylodonor" in manifest["tool_versions"]


def test_manifest_accumulates_across_steps(tmp_path):
    """Re-running one step must not erase what the others recorded."""
    cfg = make_config(tmp_path, minimize=False)
    out_dir = holo_run.build(cfg, allow_download=False)["out_dir"]
    assert set(json.loads((out_dir / "manifest.json").read_text())) >= {"transfer", "model_id"}

    holo_run.annotate_contacts(cfg)
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert "transfer" in manifest and "contacts" in manifest

    holo_run.build(cfg, allow_download=False)
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert "contacts" in manifest, "rebuilding dropped the contact record"


def test_contacts_require_a_built_model(tmp_path):
    cfg = make_config(tmp_path, minimize=False)
    with pytest.raises(ValueError, match="run `holo build` first"):
        holo_run.annotate_contacts(cfg)


def test_missing_target_structure_is_explained(tmp_path):
    path = tmp_path / "cfg.yaml"
    path.write_text(
        textwrap.dedent("""
        project: {name: demo}
        holo:
          cofactors: [{resname: LIG, donor_pdb: 1ABC}]
        """)
    )
    cfg = config_mod.load(path, repo_root=tmp_path)
    with pytest.raises(ValueError, match="no input structure"):
        holo_run.build(cfg, allow_download=False)


def test_unknown_step_is_rejected(tmp_path):
    cfg = make_config(tmp_path, minimize=False)
    with pytest.raises(ValueError, match="unknown holo step"):
        holo_run.run_step(cfg, "polish")


def test_generous_cutoff_finds_at_least_as_many_contacts(tmp_path):
    """The documented advice -- widen the cutoff instead of minimising."""
    tight = make_config(tmp_path / "a", minimize=False, cutoff=4.0)
    holo_run.build(tight, allow_download=False)
    n_tight = len(holo_run.annotate_contacts(tight)["rows"])

    wide = make_config(tmp_path / "b", minimize=False, cutoff=5.5)
    holo_run.build(wide, allow_download=False)
    n_wide = len(holo_run.annotate_contacts(wide)["rows"])

    assert n_wide >= n_tight


def test_minimisation_is_used_for_contacts_when_enabled(tmp_path):
    pytest.importorskip("openmm", reason="needs the 'dynamics' extra")
    cfg = make_config(tmp_path, minimize=True)
    out_dir = holo_run.build(cfg, allow_download=False)["out_dir"]

    result = holo_run.minimize_model(cfg)
    assert result["minimized"] is not None
    assert (out_dir / "receptor_minimized.pdb").exists()

    annotated = holo_run.annotate_contacts(cfg)
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert manifest["contacts"]["minimized"] is True
    assert manifest["contacts"]["model_used"].endswith("receptor_minimized.pdb")
    assert annotated["rows"]
