"""Enumeration of candidate grafts, under structural and positional constraints.

The job here is to be generous but not reckless. Generous, because the whole
point of the language-model stage downstream is that it can afford to triage
thousands of candidates; a stingy enumeration throws away designs before
anything has looked at them. Not reckless, because an unconstrained
enumeration will happily propose grafts straight through the active site, and
those designs waste the entire remaining budget.

Constraints, in order of importance:

1. ``exclude_positions`` and cofactor contacts -- hard no-go residues.
2. Junction placement -- seams belong in solvent-exposed loops where the
   receptor and donor backbones already superpose closely. A seam in the
   middle of a helix is a design that will not fold.
3. Length bookkeeping -- large insertions or deletions relative to the
   receptor segment change the register of everything downstream.

Correspondence between receptor and donor residues comes from the structural
superposition produced by the ``donors`` stage, not from sequence alignment.
At the identities typical of a thermophilic donor, sequence alignment is not
reliable enough to define a graft boundary.

TODO: implement. Signatures and contract are fixed.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

from phylodonor.config import Config
from phylodonor.grafting.models import Design, GraftSpec, Segment

__all__ = ["enumerate_grafts", "build_designs", "load_graft_table"]


def enumerate_grafts(
    cfg: Config,
    *,
    residue_map: Mapping[int, int],
    donor: str,
    secondary_structure: Mapping[int, str],
    relative_sasa: Mapping[int, float],
    backbone_rmsd: Mapping[int, float],
    excluded: Iterable[int] = (),
) -> list[GraftSpec]:
    """Propose single-segment grafts from one donor.

    Parameters
    ----------
    residue_map
        Receptor position -> donor position, from the structural superposition.
        Positions absent from the map are not superposable and cannot anchor a
        junction.
    secondary_structure
        Receptor position -> DSSP code. Junctions are restricted to coil/turn
        when ``junction_rules.require_loop`` is set.
    relative_sasa
        Receptor position -> relative solvent accessibility, 0-1.
    backbone_rmsd
        Receptor position -> local receptor/donor backbone RMSD in angstrom.
    excluded
        Positions that must not be modified, already unioned from
        ``candidates.exclude_positions`` and any cofactor contacts.

    Returns
    -------
    One :class:`GraftSpec` per acceptable receptor/donor segment pair, each
    carrying a short rationale.
    """
    raise NotImplementedError


def build_designs(
    grafts: Sequence[GraftSpec],
    receptor_sequence: str,
    cfg: Config,
) -> list[Design]:
    """Apply grafts to the receptor, singly and in allowed combinations.

    Combinations are limited by ``candidates.max_combinations`` and are only
    formed from non-overlapping, non-adjacent receptor segments -- two grafts
    that share a junction interact in ways neither was scored for.

    Applies grafts from the C-terminus backwards so that earlier coordinates
    stay valid as lengths change, and records the resulting junction positions
    on each :class:`Design`.
    """
    raise NotImplementedError


def load_graft_table(path: str, *, receptor: str) -> list[GraftSpec]:
    """Read hand-curated grafts from a TSV.

    Columns: ``receptor_start``, ``receptor_end``, ``donor``, ``donor_start``,
    ``donor_end``, ``rationale``.

    An escape hatch, and an important one: structural intuition and
    literature precedent routinely beat enumeration, and a curated design
    should be scored by the same machinery as a generated one so the two are
    comparable. Curated tables are project data -- keep them with your target
    config, not in the package.
    """
    raise NotImplementedError


def _is_valid_junction(position: int, cfg: Config, **evidence: Mapping[int, float] | Mapping[int, str]) -> tuple[bool, str]:
    """Whether a seam may be placed before ``position``; reason if not."""
    raise NotImplementedError
