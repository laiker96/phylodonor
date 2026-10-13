"""Sequence-guided CA superposition with iterative outlier rejection.

Why not a pure structural aligner: the donor is a crystal structure of a
different protein, often with gaps, extra domains and a different chain id.
Aligning the sequences first gives an unambiguous residue correspondence, and
superposing only the CA atoms of corresponding residues then gives a
transform that does not depend on the parts that differ.

Why outlier rejection: a global fit over every aligned pair is dominated by
the regions that genuinely differ -- a donor's extra loop drags the whole
transform and the cofactor lands in the wrong place. Rejecting pairs beyond a
cutoff and refitting converges on the conserved core, which is the frame the
cofactor site actually lives in.

Why the pocket-restricted variant: for a cofactor whose binding site is a
local pocket rather than the protein core, the globally best transform is not
the one that places that pocket best. Refitting on the donor's pocket
residues alone puts the cofactor where the receptor's equivalent pocket is,
at the cost of a worse global RMSD -- which is the right trade when the
output is cofactor coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import Superimposer

from phylodonor.holo.structure import ca_atoms, is_protein_residue, sequence_from_ca

__all__ = ["Fit", "ResidueKey", "aligned_ca_pairs", "pocket_residues", "superpose_by_sequence"]

#: ``(chain_id, residue_number, insertion_code)``; the insertion code is
#: stripped, so ``""`` means none.
ResidueKey = tuple[str, int, str]


@dataclass(frozen=True)
class Fit:
    """Outcome of a superposition. ``rmsd_a`` is over the retained pairs only."""

    rmsd_a: float
    n_pairs: int
    n_pairs_considered: int

    def as_dict(self, prefix: str) -> dict[str, float | int]:
        """Flatten into manifest fields, e.g. ``fad_rmsd_a``."""
        return {
            f"{prefix}_rmsd_a": round(self.rmsd_a, 3),
            f"{prefix}_pairs": self.n_pairs,
            f"{prefix}_pairs_considered": self.n_pairs_considered,
        }


def aligned_ca_pairs(fixed_cas: list, moving_cas: list) -> list[tuple[int, int]]:
    """Index pairs into the two CA lists, from a global sequence alignment.

    Scoring is a plain match/mismatch scheme with a steep gap-open penalty:
    the goal is a conservative one-to-one correspondence over the conserved
    core, not a biologically optimal alignment. Gapped stretches are simply
    not paired.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(sequence_from_ca(fixed_cas), sequence_from_ca(moving_cas))[0]
    pairs: list[tuple[int, int]] = []
    # The two aligned-block lists are paired by construction, so a length
    # mismatch would mean Biopython changed its contract: assert it.
    for (f0, f1), (m0, m1) in zip(alignment.aligned[0], alignment.aligned[1], strict=True):
        for offset in range(min(f1 - f0, m1 - m0)):
            pairs.append((f0 + offset, m0 + offset))
    return pairs


def _residue_key(residue) -> ResidueKey:
    return (
        residue.get_parent().id,
        int(residue.id[1]),
        str(residue.id[2]).strip(),
    )


def pocket_residues(
    structure,
    protein_chain: str,
    ligand_resname: str,
    cutoff_a: float,
    *,
    ligand_chain: str | None = None,
) -> set[ResidueKey]:
    """Protein residues with any heavy atom within ``cutoff_a`` of a ligand.

    Hydrogens are ignored on both sides: crystal structures frequently lack
    them, so including them would make the pocket depend on whether the donor
    was refined with riding hydrogens.
    """
    lig_coords = [
        atom.get_coord()
        for chain in structure[0]
        if ligand_chain is None or chain.id == ligand_chain
        for res in chain
        if res.get_resname().strip() == ligand_resname
        for atom in res
        if (atom.element or "").strip().upper() != "H"
    ]
    if not lig_coords:
        raise ValueError(f"no heavy atoms found for ligand {ligand_resname!r}")
    ligand = np.asarray(lig_coords)

    pocket: set[ResidueKey] = set()
    for res in structure[0][protein_chain]:
        if not is_protein_residue(res):
            continue
        coords = np.array(
            [a.get_coord() for a in res if (a.element or "").strip().upper() != "H"],
            dtype=float,
        )
        if coords.size == 0:
            continue
        if (
            float(np.min(np.linalg.norm(coords[:, None, :] - ligand[None, :, :], axis=2)))
            <= cutoff_a
        ):
            pocket.add(_residue_key(res))
    return pocket


def superpose_by_sequence(
    fixed_structure,
    moving_structure,
    fixed_chain: str,
    moving_chain: str,
    *,
    moving_residue_filter: set[ResidueKey] | None = None,
    min_pairs: int = 12,
    cycles: int = 8,
    outlier_cutoff_a: float = 2.5,
) -> Fit:
    """Superpose ``moving_structure`` onto ``fixed_structure``, in place.

    Every atom of ``moving_structure`` is transformed, not just the CAs used
    for fitting -- that is how the cofactor comes along.

    Parameters
    ----------
    moving_residue_filter
        Restrict fitting to these donor residues, e.g. a binding pocket from
        :func:`pocket_residues`. Pass ``cycles=1`` and a large
        ``outlier_cutoff_a`` with it: a pocket is small enough that outlier
        rejection would eat most of it.
    min_pairs
        Refuse to fit on fewer pairs than this. A transform from a handful of
        CAs is arbitrary, and silently returning it would put a cofactor
        somewhere plausible-looking and wrong.

    Raises
    ------
    ValueError
        If fewer than ``min_pairs`` usable pairs survive filtering.
    """
    fixed_cas = ca_atoms(fixed_structure, fixed_chain)
    moving_cas = ca_atoms(moving_structure, moving_chain)
    if not fixed_cas or not moving_cas:
        raise ValueError(
            f"no CA atoms to fit (fixed chain {fixed_chain!r}: {len(fixed_cas)}, "
            f"moving chain {moving_chain!r}: {len(moving_cas)})"
        )

    fixed_atoms, moving_atoms = [], []
    for fi, mi in aligned_ca_pairs(fixed_cas, moving_cas):
        if fi >= len(fixed_cas) or mi >= len(moving_cas):
            continue
        if (
            moving_residue_filter is not None
            and _residue_key(moving_cas[mi].get_parent()) not in moving_residue_filter
        ):
            continue
        fixed_atoms.append(fixed_cas[fi])
        moving_atoms.append(moving_cas[mi])

    considered = len(fixed_atoms)
    if considered < min_pairs:
        detail = " after pocket filtering" if moving_residue_filter is not None else ""
        raise ValueError(
            f"too few aligned CA pairs for superposition{detail}: {considered} < {min_pairs}"
        )

    sup = Superimposer()
    keep = list(range(considered))
    for _ in range(max(cycles, 1)):
        sup.set_atoms([fixed_atoms[i] for i in keep], [moving_atoms[i] for i in keep])
        sup.apply(moving_structure.get_atoms())
        refreshed = [
            i for i in range(considered) if (fixed_atoms[i] - moving_atoms[i]) <= outlier_cutoff_a
        ]
        # Stop when the retained set is stable, or when shrinking it further
        # would drop below the floor -- in which case the current fit stands.
        if len(refreshed) == len(keep) or len(refreshed) < min_pairs:
            break
        keep = refreshed

    sup.set_atoms([fixed_atoms[i] for i in keep], [moving_atoms[i] for i in keep])
    sup.apply(moving_structure.get_atoms())
    return Fit(rmsd_a=float(sup.rms), n_pairs=len(keep), n_pairs_considered=considered)
