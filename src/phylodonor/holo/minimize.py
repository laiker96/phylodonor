"""Restrained minimisation of a transferred holo model.

What this is for, and what it is not
------------------------------------
Transferring a cofactor from a different protein leaves small steric
overlaps where the receptor's side chains do not quite match the donor's.
This relieves them. It is a local, deterministic descent to a nearby energy
minimum -- seconds to minutes on a CPU, no trajectory, no sampling.

It is **not** molecular dynamics and it is not a stability calculation. The
energies it reports are useful only for comparing the same structure before
and after relaxation; they say nothing about whether a design is stable. For
why simulation is not on the critical path at all, see ``docs/roadmap.md``.

Why the cofactor is held fixed
------------------------------
General force fields have no parameters for an arbitrary cofactor, and
generating them (GAFF/OpenFF, AM1-BCC charges) is a slow, heavyweight step
that would dominate this stage's cost and dependencies.

It is also unnecessary here. The cofactor's internal geometry comes from a
crystal structure and does not need refining -- what needs fixing is the
receptor's local clashes *with* it. So the cofactor atoms enter the system as
fixed particles carrying element-based Lennard-Jones parameters and no
charge: enough to push overlapping protein atoms away, with nothing to
parameterise.

The honest consequence: this relieves **steric** strain only. There is no
electrostatic steering toward the cofactor, so do not read the result as an
optimised binding geometry. It is the transferred geometry with the clashes
taken out.

Requires the ``dynamics`` extra (``pip install -e ".[dynamics]"``). The import
is deferred so that the rest of the stage works without OpenMM installed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from phylodonor.holo.structure import AtomRecord, parse_atom_records, write_atom_line

__all__ = ["MinimizationResult", "MinimizationUnavailable", "minimize_holo"]


class MinimizationUnavailable(RuntimeError):
    """OpenMM is not installed, or the model cannot be parameterised."""


#: Element-based Lennard-Jones parameters for the fixed cofactor particles,
#: as (sigma_nm, epsilon_kJ_per_mol). Amber-like values; they set an
#: excluded volume per element, which is all a frozen body needs.
_LJ_BY_ELEMENT: dict[str, tuple[float, float]] = {
    "H": (0.2470, 0.0657),
    "C": (0.3400, 0.3598),
    "N": (0.3250, 0.7113),
    "O": (0.2960, 0.8786),
    "P": (0.3740, 0.8368),
    "S": (0.3564, 1.0460),
    "F": (0.3118, 0.2552),
    "CL": (0.3470, 1.1087),
    "SE": (0.3730, 1.2552),
}
_LJ_DEFAULT = _LJ_BY_ELEMENT["C"]

#: kcal/mol/A^2 -> kJ/mol/nm^2.
_KCAL_A2_TO_KJ_NM2 = 418.4


@dataclass
class MinimizationResult:
    """Outcome of one minimisation."""

    minimized_pdb: Path
    initial_energy_kj: float
    final_energy_kj: float
    max_ca_displacement_a: float
    n_protein_atoms: int
    n_fixed_ligand_atoms: int

    def as_manifest(self) -> dict[str, object]:
        return {
            "minimized_pdb": str(self.minimized_pdb),
            "initial_energy_kj_per_mol": round(self.initial_energy_kj, 2),
            "final_energy_kj_per_mol": round(self.final_energy_kj, 2),
            "max_ca_displacement_a": round(self.max_ca_displacement_a, 3),
            "n_protein_atoms": self.n_protein_atoms,
            "n_fixed_ligand_atoms": self.n_fixed_ligand_atoms,
            "note": (
                "Steric relaxation with the cofactor held fixed and uncharged. "
                "Energies are comparable only before/after this run."
            ),
        }


def _require_openmm():
    try:
        import openmm
        from openmm import app, unit  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise MinimizationUnavailable(
            "minimisation needs OpenMM. Install the extra:\n"
            '    pip install -e ".[dynamics]"\n'
            "Or disable it with holo.minimize.enabled: false -- the contact table "
            "does not require a minimised model (see docs/roadmap.md)."
        ) from exc
    return openmm, app, unit


def _split_records(records: list[AtomRecord], protein_chain: str):
    """Partition into the protein chain to relax and everything to freeze."""

    def is_protein(rec: AtomRecord) -> bool:
        return rec.kind == "ATOM" and rec.chain == protein_chain

    return [r for r in records if is_protein(r)], [r for r in records if not is_protein(r)]


def _with_terminal_oxt(protein: list[AtomRecord]) -> list[AtomRecord]:
    """Add the second carboxylate oxygen to the C-terminus if it is absent.

    Force fields match the terminal residue against a charged-terminus
    template that requires OXT. Most predicted structures include it, but a
    chain extracted from the middle of a model, or a truncated one, will not
    -- and the resulting template error is opaque relative to its cause.

    Placed in the CA-C-O plane at roughly 120 degrees from both, which is
    where a carboxylate oxygen belongs. Minimisation refines it; the point is
    only to give the template something to match.
    """
    import numpy as np

    last_resi = protein[-1].resi
    tail = [r for r in protein if r.resi == last_resi]
    if any(r.name in {"OXT", "OT2", "O2"} for r in tail):
        return protein

    by_name = {r.name: r for r in tail}
    if not {"C", "O", "CA"} <= by_name.keys():
        # Nothing sensible to build from; let the force field report it.
        return protein

    c, o, ca = by_name["C"].coord, by_name["O"].coord, by_name["CA"].coord

    def unit(v):
        norm = float(np.linalg.norm(v))
        return v / norm if norm else v

    direction = unit(-(unit(ca - c) + unit(o - c)))
    oxt = AtomRecord(
        kind="ATOM",
        name="OXT",
        resn=by_name["C"].resn,
        chain=by_name["C"].chain,
        resi=last_resi,
        icode=by_name["C"].icode,
        elem="O",
        coord=c + direction * 1.25,
        occ=by_name["C"].occ,
        bfac=by_name["C"].bfac,
    )
    return [*protein, oxt]


def minimize_holo(
    holo_pdb: Path,
    out_pdb: Path,
    *,
    protein_chain: str = "A",
    forcefield_files: tuple[str, ...] = ("amber14-all.xml",),
    max_iterations: int = 5000,
    restrain_backbone: bool = True,
    restraint_k_kcal_per_a2: float = 100.0,
    tolerance_kj_per_mol: float = 10.0,
) -> MinimizationResult:
    """Relax the protein of ``holo_pdb`` around its fixed cofactors.

    Parameters
    ----------
    restrain_backbone
        Harmonically restrain CA atoms to their input positions. On by
        default: the input fold is a prediction that has already been
        validated, and letting it drift during a clash-relief step would
        change the thing every downstream coordinate refers to.
    restraint_k_kcal_per_a2
        Restraint stiffness in kcal/mol/A^2, the unit restrained minimisation
        is conventionally specified in. Converted internally to OpenMM's
        kJ/mol/nm^2, which is a factor of 418.4 larger -- getting that
        conversion wrong gives a restraint tens of times too soft, and the
        only symptom is a backbone that quietly drifts.

    Raises
    ------
    MinimizationUnavailable
        If OpenMM is missing, or the protein cannot be parameterised -- most
        often a non-standard residue or a missing side-chain atom, which the
        error message names.
    """
    openmm, app, unit = _require_openmm()

    records = parse_atom_records(holo_pdb)
    if any("implicit/" in f for f in forcefield_files):
        raise MinimizationUnavailable(
            "implicit solvent is not supported here. The cofactors are appended to the "
            "System as fixed particles after it is built, and an implicit-solvent force "
            "is constructed with the original particle count, so the two cannot be "
            "combined without inventing GB parameters for every cofactor atom.\n"
            "Minimise in vacuum instead (the default). Keep max_iterations modest: with "
            "the backbone restrained this relieves the clashes it is meant to, but run to "
            "full convergence a vacuum minimisation over-attracts surface charges."
        )
    protein_records, ligand_records = _split_records(records, protein_chain)
    if not protein_records:
        raise MinimizationUnavailable(
            f"{holo_pdb} has no ATOM records on chain {protein_chain!r} to minimise"
        )

    protein_records = _with_terminal_oxt(protein_records)

    # Parameterise the protein alone; the cofactors are appended afterwards as
    # fixed particles, so they never reach the force field's templates.
    protein_only = out_pdb.parent / f"{out_pdb.stem}_protein_in.pdb"
    protein_only.parent.mkdir(parents=True, exist_ok=True)
    protein_only.write_text(
        "\n".join(
            [write_atom_line(i, rec) for i, rec in enumerate(protein_records, start=1)]
            + ["TER", "END"]
        )
        + "\n"
    )

    pdb = app.PDBFile(str(protein_only))
    forcefield = app.ForceField(*forcefield_files)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    try:
        modeller.addHydrogens(forcefield)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            rigidWater=False,
        )
    except (ValueError, KeyError) as exc:
        raise MinimizationUnavailable(
            f"could not parameterise the protein in {holo_pdb} with "
            f"{', '.join(forcefield_files)}: {exc}. "
            "Non-standard residues and missing side-chain atoms are the usual causes."
        ) from exc

    positions = list(modeller.positions)
    n_protein_particles = system.getNumParticles()

    nonbonded = next(
        (
            system.getForce(i)
            for i in range(system.getNumForces())
            if isinstance(system.getForce(i), openmm.NonbondedForce)
        ),
        None,
    )
    if nonbonded is None:  # pragma: no cover - every protein force field has one
        raise MinimizationUnavailable("force field produced no NonbondedForce")

    # Fixed particles: mass 0 pins them, so the protein relaxes around a rigid
    # cofactor whose crystal geometry is left exactly as transferred.
    for rec in ligand_records:
        sigma, epsilon = _LJ_BY_ELEMENT.get(rec.elem.upper(), _LJ_DEFAULT)
        system.addParticle(0.0 * unit.dalton)
        nonbonded.addParticle(0.0, sigma * unit.nanometer, epsilon * unit.kilojoule_per_mole)
        positions.append(openmm.Vec3(*(rec.coord / 10.0)) * unit.nanometer)

    if restrain_backbone:
        restraint = openmm.CustomExternalForce("0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
        # 1 kcal/mol/A^2 = 4.184 kJ/mol / (0.1 nm)^2 = 418.4 kJ/mol/nm^2
        k_kj_per_nm2 = restraint_k_kcal_per_a2 * _KCAL_A2_TO_KJ_NM2
        restraint.addGlobalParameter(
            "k", k_kj_per_nm2 * unit.kilojoule_per_mole / unit.nanometer**2
        )
        for name in ("x0", "y0", "z0"):
            restraint.addPerParticleParameter(name)
        for atom in modeller.topology.atoms():
            if atom.name == "CA":
                pos = positions[atom.index]
                restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
        system.addForce(restraint)

    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)

    initial = context.getState(getEnergy=True).getPotentialEnergy()
    openmm.LocalEnergyMinimizer.minimize(
        context,
        tolerance_kj_per_mol * unit.kilojoule_per_mole / unit.nanometer,
        max_iterations,
    )
    final_state = context.getState(getPositions=True, getEnergy=True)
    final = final_state.getPotentialEnergy()
    moved = final_state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

    # Map minimised coordinates back onto the input records. Hydrogens added
    # by the modeller are dropped: the inputs had none, and the downstream
    # consumers re-add them under their own force field.
    by_index = {}
    for atom in modeller.topology.atoms():
        if atom.element is not None and atom.element.symbol == "H":
            continue
        by_index[atom.index] = atom

    max_ca_shift = 0.0
    out_records: list[AtomRecord] = []
    heavy_indices = sorted(by_index)
    if len(heavy_indices) != len(protein_records):
        raise MinimizationUnavailable(
            f"atom-count mismatch after minimisation: {len(heavy_indices)} heavy atoms "
            f"from the force field vs {len(protein_records)} input records. "
            "This usually means the input had missing or duplicated atoms."
        )
    for index, rec in zip(heavy_indices, protein_records, strict=True):
        new = moved[index]
        if rec.name == "CA":
            max_ca_shift = max(max_ca_shift, float(((new - rec.coord) ** 2).sum() ** 0.5))
        out_records.append(
            AtomRecord(
                kind=rec.kind,
                name=rec.name,
                resn=rec.resn,
                chain=rec.chain,
                resi=rec.resi,
                icode=rec.icode,
                elem=rec.elem,
                coord=new,
                occ=rec.occ,
                bfac=rec.bfac,
            )
        )

    lines: list[str] = []
    serial = 1
    for rec in out_records:
        lines.append(write_atom_line(serial, rec))
        serial += 1
    lines.append("TER")
    for rec in ligand_records:
        lines.append(write_atom_line(serial, rec))
        serial += 1
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n")
    protein_only.unlink(missing_ok=True)

    return MinimizationResult(
        minimized_pdb=out_pdb,
        initial_energy_kj=initial.value_in_unit(unit.kilojoule_per_mole),
        final_energy_kj=final.value_in_unit(unit.kilojoule_per_mole),
        max_ca_displacement_a=max_ca_shift,
        n_protein_atoms=n_protein_particles,
        n_fixed_ligand_atoms=len(ligand_records),
    )
