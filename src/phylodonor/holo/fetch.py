"""Retrieval of donor crystal structures from the RCSB PDB.

Isolated in its own module so that everything else in the stage is pure
computation over local files, and so tests can substitute a local fixture
without monkey-patching the transfer logic.

Downloads are cached per run and never overwritten: a donor crystal is
immutable, and silently re-fetching one would make a rerun non-reproducible
if the entry were ever revised.
"""

from __future__ import annotations

import urllib.error
import urllib.request
from pathlib import Path

__all__ = ["DonorUnavailable", "download_pdb", "resolve_donor"]

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
_TIMEOUT_S = 60


class DonorUnavailable(RuntimeError):
    """A donor crystal could not be found locally or fetched."""


def download_pdb(pdb_id: str, dst: Path, *, timeout_s: int = _TIMEOUT_S) -> Path:
    """Fetch one PDB entry to ``dst``, atomically.

    Written to a temporary file and renamed, so an interrupted download
    cannot leave a truncated structure in the cache where the next run would
    treat it as valid.
    """
    url = RCSB_URL.format(pdb_id=pdb_id)
    try:
        with urllib.request.urlopen(url, timeout=timeout_s) as response:
            data = response.read()
    except urllib.error.URLError as exc:
        raise DonorUnavailable(f"could not fetch {pdb_id} from {url}: {exc}") from exc

    # An error page or a redirect to one is still a 200 with a body; checking
    # for the PDB header is cheaper and more reliable than trusting the status.
    if not data.startswith(b"HEADER"):
        raise DonorUnavailable(
            f"{url} did not return PDB-format data (no HEADER record). "
            f"Check that {pdb_id!r} is a valid PDB id."
        )

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(dst)
    return dst


def resolve_donor(
    pdb_id: str,
    cache_dir: Path,
    *,
    explicit: Path | None = None,
    allow_download: bool = True,
) -> tuple[Path, str]:
    """Locate a donor crystal, returning ``(path, provenance)``.

    Resolution order: an explicit path, then the run's cache, then RCSB. The
    provenance string goes into the stage manifest, because "where did this
    donor come from" is the first question asked of a suspicious model.
    """
    if explicit is not None:
        explicit = Path(explicit)
        if not explicit.exists():
            raise DonorUnavailable(f"explicit donor template does not exist: {explicit}")
        return explicit.resolve(), "explicit"

    cached = cache_dir / f"{pdb_id}.pdb"
    if cached.exists():
        return cached, "cached"

    if allow_download:
        return download_pdb(pdb_id, cached), "rcsb"

    raise DonorUnavailable(
        f"donor crystal {pdb_id} is not cached at {cached} and downloads are disabled. "
        f"Provide it explicitly, pre-stage it there, or enable downloads."
    )
