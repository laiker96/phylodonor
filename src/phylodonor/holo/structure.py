"""Structure parsing, normalisation and PDB record I/O.

Deliberately low-level. Biopython is used for parsing and superposition, but
the final models are written by assembling PDB lines directly, because the
downstream consumers (Amber's ``tleap``, contact-calling tools) are strict
about column positions, chain assignment and atom naming in ways that a
generic writer does not reliably satisfy.

Coordinates here are whatever the input file uses. Residue numbering is
1-based as PDB files number it, and is never renumbered -- downstream stages
map receptor positions through these files, so renumbering would silently
invalidate every graft coordinate.
"""

from __future__ import annotations

import contextlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa

__all__ = [
    "AtomRecord",
    "AMBER_TO_STANDARD",
    "STANDARD_PROTEIN_RESNAMES",
    "SOLVENT_AND_IONS",
    "VDW_RADII",
    "atom_element",
    "is_protein_residue",
    "ca_atoms",
    "model_id_from_path",
    "normalise_protein_chain",
    "parse_atom_records",
    "parse_structure",
    "first_protein_chain",
    "sequence_from_ca",
    "standard_resname",
    "write_atom_line",
    "write_model",
    "write_standard_resnames",
    "write_hydrogen_stripped",
]

# The 20 standard residues plus the protonation-state variants that Amber
# emits. Both sets are force-field vocabulary, not project data.
STANDARD_PROTEIN_RESNAMES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ASH",
        "CYM",
        "CYX",
        "GLH",
        "HID",
        "HIE",
        "HIP",
        "LYN",
    }
)

#: Amber protonation variants mapped back to their standard names, for tools
#: that only understand the canonical 20.
AMBER_TO_STANDARD = {
    "ASH": "ASP",
    "CYM": "CYS",
    "CYX": "CYS",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "LYN": "LYS",
}

SOLVENT_AND_IONS = frozenset(
    {
        "HOH",
        "WAT",
        "DOD",
        "NA",
        "CL",
        "K",
        "MG",
        "CA",
        "ZN",
        "Na+",
        "Cl-",
        "K+",
        "Mg2+",
    }
)

#: Bondi van der Waals radii in angstrom, used for clash detection.
VDW_RADII = {
    "H": 1.10,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "P": 1.80,
    "S": 1.80,
    "SE": 1.90,
    "F": 1.47,
    "CL": 1.75,
}

_STRUCTURE_SUFFIXES = (".pdb", ".ent", ".cif", ".mmcif")


@dataclass
class AtomRecord:
    """One parsed PDB ATOM/HETATM line.

    Mutable on purpose: transfer applies coordinate transforms in place.
    """

    kind: str
    name: str
    resn: str
    chain: str
    resi: str
    icode: str
    elem: str
    coord: np.ndarray
    occ: float = 1.0
    bfac: float = 20.0

    @property
    def is_hydrogen(self) -> bool:
        return self.elem in {"H", "D"}


def standard_resname(resname: str) -> str:
    """Map an Amber protonation variant back to its standard residue name."""
    cleaned = resname.strip().upper()
    return AMBER_TO_STANDARD.get(cleaned, cleaned)


def is_protein_residue(residue) -> bool:
    """Whether a Biopython residue is an amino acid, Amber variants included."""
    return (
        is_aa(residue, standard=True)
        or standard_resname(residue.get_resname()) in STANDARD_PROTEIN_RESNAMES
    )


def atom_element(line: str) -> str:
    """Element symbol for a PDB line, falling back to the atom name.

    The element columns (77-78) are optional in practice and frequently blank
    in predicted-structure output, so the fallback matters.
    """
    if len(line) >= 78:
        elem = line[76:78].strip().upper()
        if elem:
            return elem
    letters = "".join(c for c in line[12:16].strip() if c.isalpha())
    return (letters[:1] or "C").upper()


def model_id_from_path(path: Path | str) -> str:
    """A filesystem-safe identifier derived from a structure filename."""
    stem = Path(path).name
    lowered = stem.lower()
    for suffix in _STRUCTURE_SUFFIXES:
        if lowered.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    cleaned = "".join(c if c.isalnum() or c in "._-" else "_" for c in stem).strip("._-")
    return cleaned or "holo_model"


def parse_structure(path: Path):
    """Parse a PDB or mmCIF file into a Biopython structure."""
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    if suffix in {".pdb", ".ent"}:
        return PDBParser(QUIET=True).get_structure(path.stem, str(path))
    raise ValueError(
        f"unsupported structure format: {path} (expected one of {_STRUCTURE_SUFFIXES})"
    )


def first_protein_chain(structure) -> str:
    """The first chain containing any amino acid.

    Used for donor crystals, where the protein chain id varies between
    depositions and is not worth configuring per donor.
    """
    for chain in structure[0]:
        if any(is_protein_residue(res) for res in chain):
            return chain.id
    raise ValueError("no protein chain found in structure")


def ca_atoms(structure, chain_id: str) -> list:
    """CA atoms of the protein residues in one chain, in sequence order."""
    out = []
    for chain in structure[0]:
        if chain.id != chain_id:
            continue
        for res in chain:
            if is_protein_residue(res) and "CA" in res:
                out.append(res["CA"])
    return out


def sequence_from_ca(cas: list) -> str:
    """One-letter sequence for a list of CA atoms, ``X`` for anything unknown."""
    letters = []
    for ca in cas:
        resn = standard_resname(ca.get_parent().get_resname()).capitalize()
        letters.append(protein_letters_3to1.get(resn, "X"))
    return "".join(letters)


class _ProteinChainSelect(Select):
    """Keep only the standard-occupancy protein residues of one chain."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue) -> bool:
        return residue.id[0] == " " and is_protein_residue(residue)


def normalise_protein_chain(src: Path, chain_id: str, dst: Path, *, out_chain: str = "A") -> Path:
    """Extract one protein chain, canonicalise residue names, relabel the chain.

    Three normalisations, each needed by something downstream:

    - only one protein chain survives, since the rest of the pipeline works on
      a monomer and a second copy would silently double every contact count;
    - Amber protonation variants are mapped to standard names, so donor
      superposition and sequence alignment see a canonical sequence;
    - the chain is relabelled to ``out_chain`` so that cofactor chains can be
      assigned deterministically around it.
    """
    structure = parse_structure(src)
    model = next(structure.get_models())
    available = {chain.id for chain in model}
    if chain_id not in available:
        raise ValueError(
            f"chain {chain_id!r} not found in {src}; available: {', '.join(sorted(available))}"
        )

    dst.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(dst), _ProteinChainSelect(chain_id))

    lines = []
    for raw in dst.read_text().splitlines():
        if raw.startswith("ATOM"):
            resn = standard_resname(raw[17:20])
            raw = f"{raw[:17]}{resn:>3}{raw[20:21]}{out_chain}{raw[22:]}"
        elif raw.startswith("TER") and len(raw) > 21:
            raw = f"{raw[:21]}{out_chain}{raw[22:]}"
        lines.append(raw)
    dst.write_text("\n".join(lines) + "\n")
    return dst


def parse_atom_records(path: Path) -> list[AtomRecord]:
    """Read every ATOM/HETATM line of a PDB file into :class:`AtomRecord`s."""
    records = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        records.append(
            AtomRecord(
                kind=line[:6].strip(),
                name=line[12:16].strip(),
                resn=line[17:20].strip(),
                chain=line[21:22].strip(),
                resi=line[22:26].strip(),
                icode=line[26:27].strip(),
                elem=atom_element(line),
                coord=np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                ),
                occ=float(line[54:60]) if line[54:60].strip() else 1.0,
                bfac=float(line[60:66]) if line[60:66].strip() else 20.0,
            )
        )
    return records


def write_atom_line(serial: int, rec: AtomRecord) -> str:
    """Format one ATOM/HETATM line at fixed PDB column positions.

    Atom names of fewer than four characters are indented by one, as the
    format requires; four-character names start in column 13. Getting this
    wrong makes ``tleap`` misread element types.
    """
    atom_name = f"{rec.name:<4s}" if len(rec.name) >= 4 else f" {rec.name:<3s}"
    return (
        f"{rec.kind:<6}{serial:5d} {atom_name} {rec.resn:>3s} {rec.chain:1s}"
        f"{int(rec.resi):4d}{rec.icode:1s}   "
        f"{rec.coord[0]:8.3f}{rec.coord[1]:8.3f}{rec.coord[2]:8.3f}"
        f"{rec.occ:6.2f}{rec.bfac:6.2f}          {rec.elem:>2s}"
    )


def write_model(
    protein_pdb: Path,
    ligands: list[AtomRecord],
    out_pdb: Path,
    *,
    protein_chain: str = "A",
) -> Path:
    """Write a protein chain followed by ligand records, serials renumbered.

    A single ``TER`` separates protein from ligands; ligand records keep the
    chain ids assigned during transfer.
    """
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    serial = 1
    for raw in protein_pdb.read_text().splitlines():
        if raw.startswith("ATOM"):
            line = f"{raw[:6]}{serial:5d}{raw[11:]}"
            if len(line) > 21:
                line = f"{line[:21]}{protein_chain}{line[22:]}"
            lines.append(line.rstrip())
            serial += 1
        elif raw.startswith("TER"):
            lines.append("TER")
    if lines and lines[-1] != "TER":
        lines.append("TER")
    for rec in ligands:
        lines.append(write_atom_line(serial, rec))
        serial += 1
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n")
    return out_pdb


def disulfide_pairs(
    records: list[AtomRecord], cutoff_a: float = 2.3
) -> list[tuple[AtomRecord, AtomRecord]]:
    """Cysteine SG pairs within ``cutoff_a``, greedily matched nearest-first.

    Each SG is used at most once, which prevents one cysteine in a cluster
    from being reported in two bonds.
    """
    sg = [r for r in records if r.resn in {"CYS", "CYX"} and r.name == "SG"]
    pairs: list[tuple[AtomRecord, AtomRecord]] = []
    used: set[int] = set()
    for i, left in enumerate(sg):
        if i in used:
            continue
        best_j, best_d = None, cutoff_a
        for j in range(i + 1, len(sg)):
            if j in used or left.chain != sg[j].chain:
                continue
            d = float(np.linalg.norm(left.coord - sg[j].coord))
            if d <= best_d:
                best_j, best_d = j, d
        if best_j is not None:
            used.update({i, best_j})
            pairs.append((left, sg[best_j]))
    return pairs


def write_standard_resnames(src: Path, dst: Path) -> Path:
    """Rewrite Amber residue names as standard ones, with SSBOND records.

    Tools that only know the canonical 20 residues need this. The SSBOND
    header is added because mapping ``CYX`` to ``CYS`` erases the only
    indication that those cysteines were bonded.
    """
    bonds = disulfide_pairs(parse_atom_records(src))
    lines = [
        f"SSBOND {idx:3d} CYS {left.chain:1s}{int(left.resi):4d}    "
        f"CYS {right.chain:1s}{int(right.resi):4d}"
        for idx, (left, right) in enumerate(bonds, start=1)
    ]
    for raw in src.read_text().splitlines():
        if raw.startswith(("ATOM  ", "HETATM")):
            resn = raw[17:20].strip().upper()
            if resn in AMBER_TO_STANDARD:
                raw = f"{raw[:17]}{AMBER_TO_STANDARD[resn]:>3s}{raw[20:]}"
        lines.append(raw.rstrip())
    dst.write_text("\n".join(lines) + "\n")
    return dst


def write_hydrogen_stripped(src: Path, dst: Path) -> Path:
    """Write a heavy-atom-only copy, keeping CONECT records consistent.

    Dropping atoms without pruning CONECT leaves records pointing at serials
    that no longer exist, which some viewers render as spurious bonds.
    """
    removed: set[int] = set()
    kept: list[str] = []
    conect: list[str] = []
    for raw in src.read_text().splitlines():
        if raw.startswith(("ATOM  ", "HETATM")):
            if atom_element(raw) in {"H", "D"}:
                # A non-numeric serial means the line is malformed; the atom
                # is dropped either way, so there is nothing to record.
                with contextlib.suppress(ValueError):
                    removed.add(int(raw[6:11]))
                continue
            kept.append(raw.rstrip())
        elif raw.startswith("CONECT"):
            conect.append(raw)
        else:
            kept.append(raw.rstrip())

    filtered: list[str] = []
    for raw in conect:
        serials = []
        for start in range(6, len(raw), 5):
            field = raw[start : start + 5].strip()
            if field:
                # Skip unparseable CONECT fields rather than dropping the
                # whole record: the remaining serials are still valid.
                with contextlib.suppress(ValueError):
                    serials.append(int(field))
        serials = [s for s in serials if s not in removed]
        if len(serials) >= 2:
            filtered.append("CONECT" + "".join(f"{s:5d}" for s in serials))

    insert_at = len(kept)
    while insert_at > 0 and kept[insert_at - 1].startswith(("END", "MASTER")):
        insert_at -= 1
    kept[insert_at:insert_at] = filtered

    dst.write_text("\n".join(kept) + "\n")
    return dst
