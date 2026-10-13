"""Stage orchestration: transfer, optionally minimise, annotate contacts.

One entrypoint per step so that the expensive parts can be re-run
selectively, and one manifest per run recording what produced the outputs.

Minimisation is optional and off the critical path. The contact table -- the
thing ``grafting`` actually consumes -- does not need a minimised model:
relaxation shifts atoms by tenths of an angstrom, which changes distances but
very rarely changes which residues are in contact. Skipping it removes the
only heavyweight dependency in the stage. See ``docs/roadmap.md``.
"""

from __future__ import annotations

import json
import platform
import time
from pathlib import Path
from typing import Any

from phylodonor import __version__
from phylodonor.config import Config
from phylodonor.holo import contacts as contacts_mod
from phylodonor.holo import minimize as minimize_mod
from phylodonor.holo import transfer as transfer_mod
from phylodonor.holo.structure import AtomRecord, model_id_from_path, parse_atom_records

__all__ = ["build", "annotate_contacts", "minimize_model", "run_step"]

_MANIFEST_NAME = "manifest.json"


def _resolve_input(cfg: Config, input_structure: Path | None) -> Path:
    if input_structure is not None:
        return Path(input_structure).resolve()
    structure = cfg.get("target.structure", None)
    if not structure:
        raise ValueError(
            "no input structure: set target.structure in the config, or pass one "
            "explicitly. If you have no experimental structure, predict one first "
            "(see the folding helper) -- this stage relaxes and annotates a model, "
            "it does not create one."
        )
    return cfg.resolve("target.structure")


def _tool_versions() -> dict[str, str]:
    """Versions of everything that can change the numbers in the outputs."""
    versions = {"phylodonor": __version__, "python": platform.python_version()}
    for name in ("numpy", "Bio", "openmm"):
        try:
            module = __import__(name)
        except ImportError:
            continue
        versions[name] = getattr(module, "__version__", "unknown")
    return versions


def _read_manifest(out_dir: Path) -> dict[str, Any]:
    path = out_dir / _MANIFEST_NAME
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _write_manifest(out_dir: Path, updates: dict[str, Any]) -> Path:
    """Merge ``updates`` into the run manifest.

    Merged rather than replaced so that re-running one step does not discard
    what the others recorded.
    """
    manifest = _read_manifest(out_dir)
    manifest.update(updates)
    manifest["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    manifest["tool_versions"] = _tool_versions()
    path = out_dir / _MANIFEST_NAME
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return path


def _stage_dir(cfg: Config, model_id: str) -> Path:
    out = cfg.results_dir("holo") / model_id
    out.mkdir(parents=True, exist_ok=True)
    return out


def build(
    cfg: Config,
    *,
    input_structure: Path | None = None,
    allow_download: bool = True,
) -> dict[str, Any]:
    """Transfer the configured cofactors into the target structure."""
    cfg.check_stage("holo")
    source = _resolve_input(cfg, input_structure)
    model_id = model_id_from_path(source)
    out_dir = _stage_dir(cfg, model_id)

    cofactors, annotations = transfer_mod.specs_from_config(cfg)
    result = transfer_mod.transfer_cofactors(
        source,
        cofactors,
        annotations,
        out_dir,
        model_id=model_id,
        input_chain=str(cfg.get("target.input_chain", "A")),
        protein_chain=str(cfg.get("holo.protein_chain", "A")),
        allow_download=allow_download,
    )

    _write_manifest(
        out_dir,
        {
            "model_id": model_id,
            "input_structure": str(source),
            "transfer": result.as_manifest(),
        },
    )
    return {"model_id": model_id, "out_dir": out_dir, "transfer": result}


def minimize_model(cfg: Config, *, model_id: str | None = None) -> dict[str, Any]:
    """Relax the transferred model, if ``holo.minimize.enabled``."""
    cfg.check_stage("holo")
    model_id = model_id or model_id_from_path(_resolve_input(cfg, None))
    out_dir = _stage_dir(cfg, model_id)

    holo_pdb = out_dir / f"{model_id}_holo.pdb"
    if not holo_pdb.exists():
        raise ValueError(f"no transferred model at {holo_pdb}; run `holo build` first")

    if not bool(cfg.get("holo.minimize.enabled", True)):
        _write_manifest(out_dir, {"minimize": {"enabled": False}})
        return {"model_id": model_id, "out_dir": out_dir, "minimized": None}

    forcefields = cfg.get("holo.minimize.forcefield", ["amber14-all.xml"])
    if isinstance(forcefields, str):
        forcefields = [forcefields]

    result = minimize_mod.minimize_holo(
        holo_pdb,
        out_dir / f"{model_id}_minimized.pdb",
        protein_chain=str(cfg.get("holo.protein_chain", "A")),
        forcefield_files=tuple(forcefields),
        max_iterations=int(cfg.get("holo.minimize.max_iterations", 5000)),
        restrain_backbone=bool(cfg.get("holo.minimize.restrain_backbone", True)),
        restraint_k_kcal_per_a2=float(cfg.get("holo.minimize.restraint_k_kcal_per_a2", 100.0)),
    )
    _write_manifest(out_dir, {"minimize": {"enabled": True, **result.as_manifest()}})
    return {"model_id": model_id, "out_dir": out_dir, "minimized": result}


def _ligand_groups(
    model_records: list[AtomRecord],
    annotation_records: list[AtomRecord],
    cofactor_resnames: list[str],
    annotations: list[transfer_mod.AnnotationSpec],
) -> list[tuple[contacts_mod.LigandGroup, bool]]:
    """One group per cofactor, plus one per annotation-ligand set.

    The boolean says which model the group's atoms came from, and therefore
    which protein coordinates to measure against: cofactors live in the
    (possibly minimised) holo model, annotation ligands only in the
    annotation model, whose protein has not been relaxed.
    """
    groups: list[tuple[contacts_mod.LigandGroup, bool]] = []
    for resname in cofactor_resnames:
        atoms = tuple(a for a in model_records if a.resn == resname and not a.is_hydrogen)
        groups.append(
            (contacts_mod.LigandGroup(label=f"{resname.lower()}_contacts", atoms=atoms), False)
        )
    for spec in annotations:
        wanted = set(spec.resnames)
        atoms = tuple(a for a in annotation_records if a.resn in wanted and not a.is_hydrogen)
        groups.append((contacts_mod.LigandGroup(label=f"{spec.label}_contacts", atoms=atoms), True))
    return [(group, from_annotation) for group, from_annotation in groups if group.atoms]


def annotate_contacts(cfg: Config, *, model_id: str | None = None) -> dict[str, Any]:
    """Write the per-residue contact table.

    Uses the minimised model when one exists, the transferred model
    otherwise. The cutoff should be set a little generous when working from an
    unminimised model: the residue identities barely change, and for an
    exclusion list erring wide is the safe direction.
    """
    cfg.check_stage("holo")
    model_id = model_id or model_id_from_path(_resolve_input(cfg, None))
    out_dir = _stage_dir(cfg, model_id)

    minimized = out_dir / f"{model_id}_minimized.pdb"
    transferred = out_dir / f"{model_id}_holo.pdb"
    model_pdb = minimized if minimized.exists() else transferred
    if not model_pdb.exists():
        raise ValueError(f"no model at {transferred}; run `holo build` first")

    cofactors, annotations = transfer_mod.specs_from_config(cfg)
    model_records = parse_atom_records(model_pdb)

    annotation_pdb = out_dir / f"{model_id}_annotation.pdb"
    annotation_records = parse_atom_records(annotation_pdb) if annotation_pdb.exists() else []

    cutoff = float(cfg.get("holo.contacts.cutoff_angstrom", 4.0))
    protein_chain = str(cfg.get("holo.protein_chain", "A"))
    protein = [a for a in model_records if a.kind == "ATOM" and a.chain == protein_chain]

    annotation_protein = [
        a for a in annotation_records if a.kind == "ATOM" and a.chain == protein_chain
    ]

    rows: list[contacts_mod.ContactRow] = []
    groups = _ligand_groups(
        model_records, annotation_records, [c.resname for c in cofactors], annotations
    )
    for group, from_annotation in groups:
        source = annotation_protein if from_annotation else protein
        rows.extend(contacts_mod.contact_rows(source, group, cutoff, protein_chain=protein_chain))

    out_tsv = contacts_mod.write_contacts(rows, out_dir / f"{model_id}_contacts.tsv")
    counts = {
        group.label: sum(1 for r in rows if r.annotation == group.label) for group, _ in groups
    }
    severe = sorted({int(r.resi) for r in rows if r.severe_clash and r.resi.isdigit()})

    _write_manifest(
        out_dir,
        {
            "contacts": {
                "tsv": str(out_tsv),
                "model_used": str(model_pdb),
                "minimized": model_pdb == minimized,
                "cutoff_angstrom": cutoff,
                "counts": counts,
                "severe_clash_residues": severe,
            }
        },
    )
    return {
        "model_id": model_id,
        "out_dir": out_dir,
        "contacts_tsv": out_tsv,
        "rows": rows,
        "severe_clash_residues": severe,
    }


def run_step(cfg: Config, step: str, **kwargs: Any) -> dict[str, Any]:
    """Dispatch one CLI step name."""
    steps = {"build": build, "minimize": minimize_model, "qc": annotate_contacts}
    if step not in steps:
        raise ValueError(f"unknown holo step {step!r}; expected one of {', '.join(steps)}")
    return steps[step](cfg, **kwargs)
