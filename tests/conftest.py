"""Synthetic structure fixtures.

Everything here is generated, so the tests need no network, no reference
database and no real PDB entry. The geometry is a poly-alanine helix: not a
real protein, but it exercises exactly what the stage relies on -- parseable
residues, CA atoms in sequence order, and a ligand at a known position
relative to the protein.

The key trick for testing transfer: build a receptor, then build a "donor" by
applying a known rigid transform to the same coordinates with a ligand
attached. A correct superposition must undo that transform exactly, so the
transferred ligand should land back on its original coordinates. Any error in
the alignment or the coordinate handling shows up as a position error.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

# Backbone offsets from CA, roughly physical and entirely sufficient here.
_BACKBONE = {
    "N": np.array([-0.60, -1.10, -0.35]),
    "C": np.array([1.10, 0.55, 0.40]),
    "O": np.array([1.45, 0.60, 1.60]),
}
_HELIX_RADIUS_A = 2.3
_HELIX_RISE_A = 1.5
_HELIX_TURN_DEG = 100.0


def helix_ca_coords(n_residues: int) -> np.ndarray:
    """CA positions of an idealised alpha helix along z."""
    out = []
    for i in range(n_residues):
        theta = math.radians(_HELIX_TURN_DEG * i)
        out.append(
            [
                _HELIX_RADIUS_A * math.cos(theta),
                _HELIX_RADIUS_A * math.sin(theta),
                _HELIX_RISE_A * i,
            ]
        )
    return np.array(out, dtype=float)


def rotation_matrix(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotation about ``axis`` by ``angle_deg``, via Rodrigues' formula."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    theta = math.radians(angle_deg)
    kx, ky, kz = axis
    k = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]], dtype=float)
    return np.eye(3) + math.sin(theta) * k + (1 - math.cos(theta)) * (k @ k)


def _atom_line(kind, serial, name, resn, chain, resi, coord, element):
    atom_name = f"{name:<4s}" if len(name) >= 4 else f" {name:<3s}"
    return (
        f"{kind:<6}{serial:5d} {atom_name} {resn:>3s} {chain:1s}{resi:4d}    "
        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
        f"{1.00:6.2f}{20.00:6.2f}          {element:>2s}"
    )


def write_structure(
    path: Path,
    *,
    n_residues: int = 60,
    sequence: str | None = None,
    protein_chain: str = "A",
    ligands: dict[str, list[tuple[str, np.ndarray, str]]] | None = None,
    transform: tuple[np.ndarray, np.ndarray] | None = None,
) -> Path:
    """Write a synthetic PDB.

    Parameters
    ----------
    sequence
        Three-letter residue names, one per residue, cycled. Defaults to all
        alanine.
    ligands
        ``{chain_id: [(resname, coord, element), ...]}``. Each entry becomes
        one HETATM; entries sharing a resname and chain are grouped into one
        residue.
    transform
        ``(rotation, translation)`` applied to every atom, protein and ligand
        alike, so a transformed copy stays internally consistent.
    """
    resnames = ["ALA"] * n_residues
    if sequence:
        resnames = [sequence[i % len(sequence)] for i in range(n_residues)]

    cas = helix_ca_coords(n_residues)
    rotation, translation = transform if transform else (np.eye(3), np.zeros(3))

    def place(coord: np.ndarray) -> np.ndarray:
        return rotation @ np.asarray(coord, dtype=float) + translation

    lines: list[str] = []
    serial = 1
    for i, ca in enumerate(cas):
        resi = i + 1
        for name, offset in (
            ("N", _BACKBONE["N"]),
            ("CA", np.zeros(3)),
            ("C", _BACKBONE["C"]),
            ("O", _BACKBONE["O"]),
        ):
            lines.append(
                _atom_line(
                    "ATOM",
                    serial,
                    name,
                    resnames[i],
                    protein_chain,
                    resi,
                    place(ca + offset),
                    name[0],
                )
            )
            serial += 1
    lines.append("TER")

    for chain_id, atoms in (ligands or {}).items():
        grouped: dict[str, list[tuple[str, np.ndarray, str]]] = {}
        for resname, coord, element in atoms:
            grouped.setdefault(resname, []).append((resname, coord, element))
        for resi, (resname, entries) in enumerate(grouped.items(), start=1):
            for index, (_, coord, element) in enumerate(entries, start=1):
                lines.append(
                    _atom_line(
                        "HETATM",
                        serial,
                        f"{element}{index}",
                        resname,
                        chain_id,
                        resi,
                        place(coord),
                        element,
                    )
                )
                serial += 1
    lines.append("END")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def receptor_pdb(tmp_path: Path) -> Path:
    """An apo receptor: 60-residue helix, chain A, no ligands."""
    return write_structure(tmp_path / "receptor.pdb", n_residues=60)


@pytest.fixture
def ligand_coords() -> np.ndarray:
    """Four cofactor atoms sitting beside the middle of the helix."""
    centre = helix_ca_coords(60)[30]
    return np.array(
        [
            centre + [4.0, 0.0, 0.0],
            centre + [5.4, 0.0, 0.0],
            centre + [4.7, 1.3, 0.0],
            centre + [4.7, 0.0, 1.3],
        ],
        dtype=float,
    )


@pytest.fixture
def donor_transform() -> tuple[np.ndarray, np.ndarray]:
    """A non-trivial rigid transform the superposition has to undo."""
    return rotation_matrix(np.array([0.3, 0.8, 0.5]), 47.0), np.array([13.0, -7.5, 22.0])


@pytest.fixture
def donor_pdb(tmp_path, ligand_coords, donor_transform) -> Path:
    """A donor crystal: the same helix, rigidly moved, carrying a ligand."""
    return write_structure(
        tmp_path / "donor_1ABC.pdb",
        n_residues=60,
        protein_chain="X",
        ligands={"Y": [("LIG", coord, "C") for coord in ligand_coords]},
        transform=donor_transform,
    )
