"""Cofactor transfer: donor crystal coordinates into the receptor frame.

The problem this solves: structure predictors give you an apo model. For an
enzyme whose function depends on a cofactor, an apo model is the wrong input
to minimisation, to dynamics, and to any contact-based analysis -- the
cofactor pocket collapses without its occupant. Experimental holo structures
of close homologs do exist, so the cofactor's coordinates can be borrowed.

The method, per cofactor:

1. Resolve a donor crystal containing that cofactor (explicit, cached, RCSB).
2. Superpose the donor's protein onto the receptor by sequence-guided CA
   alignment with outlier rejection. Every donor atom moves, cofactor
   included.
3. Optionally refit on the donor's binding pocket alone, for a cofactor in a
   local site where the global transform is not the one that places the
   pocket best.
4. Copy the cofactor's heavy atoms, now in the receptor frame, onto a
   deterministically assigned chain.

Hydrogens are dropped: donor crystals disagree about whether they have any,
and the force field adds them during preparation.

Annotation ligands (nucleic acids, substrate analogues, inhibitors) are
handled separately. They are transferred into a second model used only for
contact annotation, never into the structure that gets minimised -- they are
evidence about where the substrate sits, not part of the enzyme.
"""

from __future__ import annotations

import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from phylodonor.holo import fetch
from phylodonor.holo.structure import (
    AtomRecord,
    first_protein_chain,
    normalise_protein_chain,
    parse_structure,
    write_model,
)
from phylodonor.holo.superpose import pocket_residues, superpose_by_sequence

__all__ = [
    "AnnotationSpec",
    "CofactorSpec",
    "TransferResult",
    "specs_from_config",
    "transfer_cofactors",
]

#: Chains for transferred ligands, assigned in config order after the
#: protein. Skips the protein's own chain so assignment stays deterministic
#: regardless of how many cofactors there are.
_LIGAND_CHAIN_POOL = string.ascii_uppercase

# Below this many CA pairs a whole-protein fit is not trustworthy. Pocket
# refits use their own, much lower floor -- a pocket is only a few residues.
_DEFAULT_MIN_FIT_PAIRS = 40
_POCKET_MIN_FIT_PAIRS = 4


@dataclass(frozen=True)
class CofactorSpec:
    """One cofactor to transfer into the receptor."""

    resname: str
    donor_pdb: str
    donor_chain: str | None = None
    """Protein chain of the donor to fit on. ``None`` picks the first one."""
    donor_ligand_chain: str | None = None
    """Chain the cofactor sits on in the donor. ``None`` searches all chains."""
    out_chain: str | None = None
    """Chain to write it on. ``None`` assigns one automatically."""
    out_resi: int = 1
    charge: int | None = None
    """Formal charge, passed through to force-field parameterisation."""
    pocket_cutoff_angstrom: float | None = None
    """Set to refit on the donor's binding pocket before copying."""
    min_fit_pairs: int = _DEFAULT_MIN_FIT_PAIRS
    explicit_path: Path | None = None


@dataclass(frozen=True)
class AnnotationSpec:
    """Ligands transferred for contact annotation only.

    ``from_donor`` names the cofactor whose already-superposed donor crystal
    these come from, so no second download or fit is needed.
    """

    label: str
    resnames: tuple[str, ...]
    from_donor: str
    preserve_chains: bool = True
    """Keep the donor's chain separation (e.g. two nucleic-acid strands),
    remapped onto fresh chain ids. ``False`` collapses them onto one."""


@dataclass
class TransferResult:
    """Outputs of a transfer run, ready to serialise into a manifest."""

    protein_pdb: Path
    holo_pdb: Path
    annotation_pdb: Path | None
    cofactor_chains: dict[str, str] = field(default_factory=dict)
    cofactor_charges: dict[str, int] = field(default_factory=dict)
    donors: dict[str, dict[str, str]] = field(default_factory=dict)
    fits: dict[str, Any] = field(default_factory=dict)
    annotation_counts: dict[str, int] = field(default_factory=dict)

    def as_manifest(self) -> dict[str, Any]:
        return {
            "protein_pdb": str(self.protein_pdb),
            "holo_pdb": str(self.holo_pdb),
            "annotation_pdb": str(self.annotation_pdb) if self.annotation_pdb else None,
            "cofactor_chains": self.cofactor_chains,
            "cofactor_charges": self.cofactor_charges,
            "donors": self.donors,
            "fits": self.fits,
            "annotation_counts": self.annotation_counts,
        }


def _as_int(value: Any, name: str) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"holo config: {name} must be an integer, got {value!r}") from exc


def specs_from_config(cfg) -> tuple[list[CofactorSpec], list[AnnotationSpec]]:
    """Build specs from the ``holo`` section of a config.

    Validates the cross-references here rather than at use time: an
    annotation ligand pointing at a cofactor that is not configured should
    fail before a donor crystal is downloaded, not after.
    """
    raw_cofactors = cfg.require("holo.cofactors")
    cofactors: list[CofactorSpec] = []
    for entry in raw_cofactors:
        if not isinstance(entry, dict) or "resname" not in entry or "donor_pdb" not in entry:
            raise ValueError(
                f"holo.cofactors entries need at least 'resname' and 'donor_pdb'; got {entry!r}"
            )
        explicit = entry.get("explicit_path")
        cofactors.append(
            CofactorSpec(
                resname=str(entry["resname"]).strip().upper(),
                donor_pdb=str(entry["donor_pdb"]).strip().upper(),
                donor_chain=entry.get("donor_chain"),
                donor_ligand_chain=entry.get("donor_ligand_chain"),
                out_chain=entry.get("out_chain"),
                out_resi=int(entry.get("out_resi", 1)),
                charge=_as_int(entry.get("charge"), "holo.cofactors[].charge"),
                pocket_cutoff_angstrom=entry.get("pocket_cutoff_angstrom"),
                min_fit_pairs=int(entry.get("min_fit_pairs", _DEFAULT_MIN_FIT_PAIRS)),
                explicit_path=Path(explicit) if explicit else None,
            )
        )

    known = {c.resname for c in cofactors}
    annotations: list[AnnotationSpec] = []
    for entry in cfg.get("holo.annotation_ligands", []) or []:
        from_donor = str(entry.get("from_donor", "")).strip().upper()
        if from_donor not in known:
            raise ValueError(
                f"holo.annotation_ligands[{entry.get('label')!r}].from_donor={from_donor!r} "
                f"is not a configured cofactor; known: {', '.join(sorted(known)) or 'none'}"
            )
        annotations.append(
            AnnotationSpec(
                label=str(entry["label"]),
                resnames=tuple(str(r).strip().upper() for r in entry["resnames"]),
                from_donor=from_donor,
                preserve_chains=bool(entry.get("preserve_chains", True)),
            )
        )
    return cofactors, annotations


def _assign_chains(specs: list[CofactorSpec], protein_chain: str) -> dict[str, str]:
    """Give every cofactor a chain id, honouring explicit requests."""
    taken = {protein_chain.upper()} | {s.out_chain.upper() for s in specs if s.out_chain}
    assigned: dict[str, str] = {}
    pool = (c for c in _LIGAND_CHAIN_POOL if c not in taken)
    for spec in specs:
        if spec.out_chain:
            assigned[spec.resname] = spec.out_chain.upper()
            continue
        try:
            assigned[spec.resname] = next(pool)
        except StopIteration as exc:  # pragma: no cover - 25 cofactors is absurd
            raise ValueError("ran out of chain ids for cofactors") from exc
    return assigned


def _extract_ligand(
    structure,
    resnames: set[str],
    *,
    out_chain: str,
    out_resi: int | None,
    ligand_chain: str | None = None,
) -> list[AtomRecord]:
    """Heavy atoms of the named residues, relabelled onto ``out_chain``.

    ``out_resi=None`` preserves the donor's residue numbering, which matters
    for multi-residue ligands such as a nucleic-acid strand.
    """
    out: list[AtomRecord] = []
    for chain in structure[0]:
        if ligand_chain is not None and chain.id != ligand_chain:
            continue
        for res in chain:
            if res.get_resname().strip().upper() not in resnames:
                continue
            for atom in res:
                elem = (atom.element or "").strip().upper()
                if elem in {"H", "D"}:
                    continue
                out.append(
                    AtomRecord(
                        kind="HETATM",
                        name=atom.get_name().strip(),
                        resn=res.get_resname().strip().upper(),
                        chain=out_chain,
                        resi=str(out_resi if out_resi is not None else res.get_id()[1]),
                        icode="" if out_resi is not None else str(res.get_id()[2]).strip(),
                        elem=elem or atom.get_name().strip()[:1].upper(),
                        coord=np.array(atom.get_coord(), dtype=float),
                        occ=atom.get_occupancy() or 1.0,
                        bfac=atom.get_bfactor() or 20.0,
                    )
                )
    return out


def _find_ligand_chain(structure, resname: str) -> str:
    for chain in structure[0]:
        for res in chain:
            if res.get_resname().strip().upper() == resname:
                return chain.id
    raise ValueError(f"donor structure does not contain ligand {resname!r}")


def _extract_annotation(
    structure,
    spec: AnnotationSpec,
    chain_pool: list[str],
) -> list[AtomRecord]:
    """Annotation ligands, each donor chain remapped onto a fresh chain id."""
    wanted = set(spec.resnames)
    donor_chains = [
        chain.id
        for chain in structure[0]
        if any(res.get_resname().strip().upper() in wanted for res in chain)
    ]
    if not donor_chains:
        return []

    out: list[AtomRecord] = []
    if spec.preserve_chains:
        for index, donor_chain in enumerate(donor_chains):
            out_chain = chain_pool[index] if index < len(chain_pool) else chain_pool[-1]
            out.extend(
                _extract_ligand(
                    structure, wanted, out_chain=out_chain, out_resi=None, ligand_chain=donor_chain
                )
            )
    else:
        out.extend(_extract_ligand(structure, wanted, out_chain=chain_pool[0], out_resi=None))
    return out


def transfer_cofactors(
    input_structure: Path,
    cofactors: list[CofactorSpec],
    annotations: list[AnnotationSpec],
    out_dir: Path,
    *,
    model_id: str,
    input_chain: str = "A",
    protein_chain: str = "A",
    allow_download: bool = True,
) -> TransferResult:
    """Build a cofactor-bound model of ``input_structure``.

    Writes ``<model_id>_holo.pdb`` (protein plus cofactors, for minimisation)
    and, when annotation ligands are configured, ``<model_id>_annotation.pdb``
    (protein plus those ligands, for contact calling only).

    Raises
    ------
    ValueError
        On an unusable configuration, a missing chain, a donor lacking the
        named cofactor, or a superposition with too few residue pairs. Every
        one of these would otherwise produce a model that looks fine and is
        wrong, so none of them are warnings.
    """
    if not cofactors:
        raise ValueError("no cofactors configured; holo.cofactors must not be empty")

    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "donor_crystals"

    protein_pdb = normalise_protein_chain(
        input_structure,
        input_chain,
        out_dir / "intermediates" / f"{model_id}_protein.pdb",
        out_chain=protein_chain,
    )
    receptor = parse_structure(protein_pdb)

    chains = _assign_chains(cofactors, protein_chain)
    used_chains = {protein_chain.upper()} | set(chains.values())
    annotation_pool = [c for c in _LIGAND_CHAIN_POOL if c not in used_chains]

    result = TransferResult(
        protein_pdb=protein_pdb,
        holo_pdb=out_dir / f"{model_id}_holo.pdb",
        annotation_pdb=None,
        cofactor_chains=chains,
    )

    cofactor_records: list[AtomRecord] = []
    # Donor structures are kept after fitting so annotation ligands can reuse
    # the same transform instead of refitting.
    fitted_donors: dict[str, Any] = {}

    for spec in cofactors:
        donor_path, provenance = fetch.resolve_donor(
            spec.donor_pdb,
            cache_dir,
            explicit=spec.explicit_path,
            allow_download=allow_download,
        )
        donor = parse_structure(donor_path)
        donor_protein_chain = spec.donor_chain or first_protein_chain(donor)

        global_fit = superpose_by_sequence(
            receptor,
            donor,
            protein_chain,
            donor_protein_chain,
            min_pairs=spec.min_fit_pairs,
        )
        result.fits.update(global_fit.as_dict(f"{spec.resname.lower()}_global"))

        if spec.pocket_cutoff_angstrom is not None:
            pocket = pocket_residues(
                donor,
                donor_protein_chain,
                spec.resname,
                float(spec.pocket_cutoff_angstrom),
                ligand_chain=spec.donor_ligand_chain,
            )
            # One cycle, no outlier rejection: the pocket is already the
            # selection, and trimming it further would leave too little to fit.
            pocket_fit = superpose_by_sequence(
                receptor,
                donor,
                protein_chain,
                donor_protein_chain,
                moving_residue_filter=pocket,
                min_pairs=_POCKET_MIN_FIT_PAIRS,
                cycles=1,
                outlier_cutoff_a=float("inf"),
            )
            result.fits.update(pocket_fit.as_dict(f"{spec.resname.lower()}_pocket"))

        ligand_chain = spec.donor_ligand_chain or _find_ligand_chain(donor, spec.resname)
        records = _extract_ligand(
            donor,
            {spec.resname},
            out_chain=chains[spec.resname],
            out_resi=spec.out_resi,
            ligand_chain=ligand_chain,
        )
        if not records:
            raise ValueError(
                f"donor {spec.donor_pdb} chain {ligand_chain!r} yielded no heavy atoms "
                f"for cofactor {spec.resname!r}"
            )
        cofactor_records.extend(records)

        if spec.charge is not None:
            result.cofactor_charges[spec.resname] = spec.charge
        result.donors[spec.resname] = {
            "pdb_id": spec.donor_pdb,
            "path": str(donor_path),
            "source": provenance,
            "protein_chain": donor_protein_chain,
            "ligand_chain": ligand_chain,
        }
        fitted_donors[spec.resname] = donor

    write_model(protein_pdb, cofactor_records, result.holo_pdb, protein_chain=protein_chain)

    if annotations:
        annotation_records: list[AtomRecord] = []
        pool = list(annotation_pool)
        for spec in annotations:
            donor = fitted_donors[spec.from_donor]
            records = _extract_annotation(donor, spec, pool)
            result.annotation_counts[spec.label] = len(records)
            if records:
                # Consume the chain ids this spec actually used, so a second
                # annotation group does not land on the same chains.
                consumed = {r.chain for r in records}
                pool = [c for c in pool if c not in consumed] or pool
            annotation_records.extend(records)

        result.annotation_pdb = out_dir / f"{model_id}_annotation.pdb"
        write_model(
            protein_pdb,
            cofactor_records + annotation_records,
            result.annotation_pdb,
            protein_chain=protein_chain,
        )

    return result
