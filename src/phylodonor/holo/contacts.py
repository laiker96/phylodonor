"""Per-residue contact annotation against cofactors and annotation ligands.

Produces the table that downstream stages depend on: which receptor residues
touch which ligand, how closely, and whether the contact is chemically
sensible or a clash. ``grafting`` reads it to build its no-go list, so the
residue numbering here is the receptor's and must stay that way.

Two judgements are encoded, and both are heuristics worth knowing about:

**van der Waals overlap** is the sum of the two Bondi radii minus the
distance. Positive means the atoms are closer than their radii allow. Small
positive overlaps are normal in real structures; large ones mean the geometry
is wrong.

**Hydrogen-bond compatibility** distinguishes the common false positive. An
N or O pair at 2.2-3.2 angstrom shows a large formal overlap while being a
perfectly good hydrogen bond. Flagging those as clashes would condemn every
correctly placed cofactor, so an overlapping pair of N/O atoms no closer than
2.2 angstrom is reported as an H-bond-compatible close contact instead.

A contact is only called severe when the overlap exceeds 0.6 angstrom *and*
the pair cannot be a hydrogen bond. That is the flag to act on.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from phylodonor.holo.structure import VDW_RADII, AtomRecord

__all__ = [
    "FIELDS",
    "ContactRow",
    "LigandGroup",
    "closest_approach",
    "contact_rows",
    "write_contacts",
]

#: Output column order. Downstream stages index by name, but a stable order
#: keeps diffs between runs readable.
FIELDS = (
    "annotation",
    "chain",
    "resi",
    "resn",
    "min_distance_a",
    "protein_atom",
    "target_chain",
    "target_resi",
    "target_resn",
    "target_atom",
    "max_vdw_overlap_a",
    "hbond_compatible",
    "severe_clash",
)

_DEFAULT_RADIUS = 1.70
_HBOND_MIN_DISTANCE_A = 2.2
_HBOND_OVERLAP_A = 0.4
_SEVERE_OVERLAP_A = 0.6
_HBOND_ELEMENTS = frozenset({"N", "O"})


@dataclass(frozen=True)
class LigandGroup:
    """A set of atoms to measure the protein against, under one label."""

    label: str
    atoms: tuple[AtomRecord, ...]


@dataclass(frozen=True)
class ContactRow:
    """One protein residue's closest approach to one ligand group."""

    annotation: str
    chain: str
    resi: str
    resn: str
    min_distance_a: float
    protein_atom: str
    target_chain: str
    target_resi: str
    target_resn: str
    target_atom: str
    max_vdw_overlap_a: float
    hbond_compatible: bool
    severe_clash: bool

    def as_row(self) -> dict[str, str]:
        return {
            "annotation": self.annotation,
            "chain": self.chain,
            "resi": self.resi,
            "resn": self.resn,
            "min_distance_a": f"{self.min_distance_a:.3f}",
            "protein_atom": self.protein_atom,
            "target_chain": self.target_chain,
            "target_resi": self.target_resi,
            "target_resn": self.target_resn,
            "target_atom": self.target_atom,
            "max_vdw_overlap_a": f"{self.max_vdw_overlap_a:.3f}",
            "hbond_compatible": "yes" if self.hbond_compatible else "no",
            "severe_clash": "yes" if self.severe_clash else "no",
        }


def _radius(element: str) -> float:
    return VDW_RADII.get(element.upper(), _DEFAULT_RADIUS)


def closest_approach(
    residue_atoms: list[AtomRecord],
    target_atoms: list[AtomRecord],
) -> dict[str, object] | None:
    """Nearest heavy-atom pair, and the worst van der Waals overlap.

    These are two different pairs in general: the closest pair is not
    necessarily the most overlapping one, because radii differ by element. The
    distance describes the contact; the overlap judges it.

    Returns ``None`` when either side has no heavy atoms.
    """
    protein = [a for a in residue_atoms if not a.is_hydrogen]
    target = [a for a in target_atoms if not a.is_hydrogen]
    if not protein or not target:
        return None

    p_coords = np.array([a.coord for a in protein], dtype=float)
    t_coords = np.array([a.coord for a in target], dtype=float)
    distances = np.linalg.norm(p_coords[:, None, :] - t_coords[None, :, :], axis=2)

    pi, ti = np.unravel_index(int(np.argmin(distances)), distances.shape)
    min_distance = float(distances[pi, ti])

    radii_sum = np.add.outer(
        np.array([_radius(a.elem) for a in protein]),
        np.array([_radius(a.elem) for a in target]),
    )
    overlaps = radii_sum - distances
    oi, oj = np.unravel_index(int(np.argmax(overlaps)), overlaps.shape)
    max_overlap = float(overlaps[oi, oj])

    worst_protein, worst_target = protein[oi], target[oj]
    worst_distance = float(distances[oi, oj])
    could_hbond = (
        worst_protein.elem in _HBOND_ELEMENTS
        and worst_target.elem in _HBOND_ELEMENTS
        and worst_distance >= _HBOND_MIN_DISTANCE_A
    )

    return {
        "min_distance_a": min_distance,
        "protein_atom": protein[pi].name,
        "target": target[ti],
        "max_vdw_overlap_a": max_overlap,
        "hbond_compatible": max_overlap > _HBOND_OVERLAP_A and could_hbond,
        "severe_clash": max_overlap > _SEVERE_OVERLAP_A and not could_hbond,
    }


def _group_by_residue(
    atoms: list[AtomRecord],
    protein_chain: str,
) -> dict[tuple[str, str, str], list[AtomRecord]]:
    grouped: dict[tuple[str, str, str], list[AtomRecord]] = {}
    for atom in atoms:
        if atom.kind != "ATOM" or atom.chain != protein_chain or atom.is_hydrogen:
            continue
        grouped.setdefault((atom.chain, atom.resi, atom.resn), []).append(atom)
    return grouped


def contact_rows(
    protein_atoms: list[AtomRecord],
    group: LigandGroup,
    cutoff_a: float,
    *,
    protein_chain: str = "A",
) -> list[ContactRow]:
    """Residues of ``protein_chain`` within ``cutoff_a`` of ``group``.

    Sorted by residue number so the output is stable across runs; residues
    with non-numeric ids sort last rather than crashing the comparison.
    """
    rows: list[ContactRow] = []
    for (chain, resi, resn), atoms in _group_by_residue(protein_atoms, protein_chain).items():
        geom = closest_approach(atoms, list(group.atoms))
        if geom is None or geom["min_distance_a"] > cutoff_a:
            continue
        target: AtomRecord = geom["target"]  # type: ignore[assignment]
        rows.append(
            ContactRow(
                annotation=group.label,
                chain=chain,
                resi=resi,
                resn=resn,
                min_distance_a=float(geom["min_distance_a"]),
                protein_atom=str(geom["protein_atom"]),
                target_chain=target.chain,
                target_resi=target.resi,
                target_resn=target.resn,
                target_atom=target.name,
                max_vdw_overlap_a=float(geom["max_vdw_overlap_a"]),
                hbond_compatible=bool(geom["hbond_compatible"]),
                severe_clash=bool(geom["severe_clash"]),
            )
        )
    rows.sort(key=lambda r: (int(r.resi) if r.resi.lstrip("-").isdigit() else 10**9, r.resn))
    return rows


def write_contacts(rows: list[ContactRow], out_tsv: Path) -> Path:
    """Write the contact table, including the header when there are no rows.

    An empty table with a header is a result ("nothing is in contact"); a
    missing file is an error. Downstream stages distinguish the two.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(row.as_row() for row in rows)
    return out_tsv
