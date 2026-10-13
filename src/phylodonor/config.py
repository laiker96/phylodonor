"""Configuration loading and validation.

Every stage takes its parameters from a single YAML file (see
``configs/template.yaml``). This module is the only place that reads it.

The design rule this enforces: no protein-specific value -- accession,
residue number, organism, melting temperature, input path -- may appear
anywhere in ``src/``. If a stage needs one, it reads it from here. That keeps
the package reusable across targets, and it keeps published code free of the
data it was originally developed against.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

__all__ = ["Config", "ConfigError", "load"]


class ConfigError(ValueError):
    """Raised when a config file is missing, malformed, or internally inconsistent."""


# Weight blocks that must sum to 1.0, as dotted paths into the config.
_WEIGHT_BLOCKS = ("donors.weights", "grafting.score.weights")

# Stages and the config sections they require to be non-null.
_STAGE_REQUIREMENTS: dict[str, tuple[str, ...]] = {
    "homologs": ("target.fasta", "search.database"),
    "thermal": ("thermal.regimes", "thermal.cluster_regimes"),
    "donors": ("donors.weights", "donors.min_ogt_delta"),
    "grafting": ("grafting.language_model.sequence_model", "grafting.score.weights"),
    "holo": ("holo.cofactors",),
    "dynamics": ("dynamics.engine", "dynamics.forcefield"),
    "scoring": ("scoring.conservation.method",),
    # Opt-in; required keys only matter if you choose to run them.
    "phylogeny": ("phylogeny.aligner",),
}

_SLUG = re.compile(r"^[a-z0-9][a-z0-9_-]*$")


@dataclass(frozen=True)
class Config:
    """A validated target configuration.

    Access nested values with :meth:`get` using dotted paths, which keeps
    stage code readable and gives a single error path for missing keys.
    """

    raw: dict[str, Any]
    path: Path
    _repo_root: Path = field(repr=False, default=Path.cwd())

    # -- access ---------------------------------------------------------------

    def get(self, dotted: str, default: Any = _SLUG) -> Any:
        """Return the value at ``dotted``, e.g. ``"donors.min_ogt_delta"``.

        Raises :class:`ConfigError` for a missing path unless ``default`` is
        given. The sentinel is deliberately odd so that ``None`` remains a
        usable default.
        """
        node: Any = self.raw
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                if default is not _SLUG:
                    return default
                raise ConfigError(f"{self.path}: missing required key '{dotted}'")
            node = node[part]
        return node

    def require(self, dotted: str) -> Any:
        """Like :meth:`get`, but also rejects ``None`` and empty collections."""
        value = self.get(dotted)
        if value is None or (hasattr(value, "__len__") and len(value) == 0):
            raise ConfigError(f"{self.path}: key '{dotted}' must be set before running this stage")
        return value

    def resolve(self, dotted: str, *, must_exist: bool = True) -> Path:
        """Resolve a config value as a path relative to the repository root."""
        value = self.require(dotted)
        resolved = (self._repo_root / str(value)).resolve()
        if must_exist and not resolved.exists():
            raise ConfigError(f"{self.path}: '{dotted}' points at {resolved}, which does not exist")
        return resolved

    # -- derived locations ----------------------------------------------------

    @property
    def name(self) -> str:
        return str(self.get("project.name"))

    def results_dir(self, stage: str) -> Path:
        """Output directory for ``stage``, created on access."""
        root = self._repo_root / str(self.get("paths.results_root", "results"))
        out = root / self.name / stage
        out.mkdir(parents=True, exist_ok=True)
        return out

    def scratch_dir(self, stage: str) -> Path:
        root = self._repo_root / str(self.get("paths.scratch", ".scratch"))
        out = root / self.name / stage
        out.mkdir(parents=True, exist_ok=True)
        return out

    # -- validation -----------------------------------------------------------

    def check_stage(self, stage: str) -> None:
        """Fail early if this config cannot drive ``stage``.

        Called at the top of each stage entrypoint so that a missing database
        path is reported in a second rather than after an hour of alignment.
        """
        if stage not in _STAGE_REQUIREMENTS:
            raise ConfigError(f"unknown stage '{stage}'")
        missing = []
        for key in _STAGE_REQUIREMENTS[stage]:
            try:
                self.require(key)
            except ConfigError:
                missing.append(key)
        if missing:
            raise ConfigError(
                f"{self.path}: stage '{stage}' needs these keys set: {', '.join(missing)}"
            )


def _check_weights(raw: dict[str, Any], path: Path) -> None:
    for block in _WEIGHT_BLOCKS:
        node: Any = raw
        for part in block.split("."):
            if not isinstance(node, dict) or part not in node:
                node = None
                break
            node = node[part]
        if not isinstance(node, dict) or not node:
            continue
        total = sum(float(v) for v in node.values())
        if abs(total - 1.0) > 1e-6:
            raise ConfigError(
                f"{path}: '{block}' must sum to 1.0, got {total:.4f}. "
                "Weights are relative, so renormalise rather than guessing."
            )


def _check_ranges(raw: dict[str, Any], path: Path) -> None:
    """Validate 1-based inclusive coordinate ranges wherever they appear."""
    for block, items in (
        ("domains.definitions", raw.get("domains", {}).get("definitions") or []),
        (
            "grafting.candidates.candidate_regions",
            (raw.get("grafting", {}).get("candidates", {}) or {}).get("candidate_regions") or [],
        ),
    ):
        for item in items:
            if not isinstance(item, dict):
                raise ConfigError(f"{path}: entries in '{block}' must be mappings")
            start, end = item.get("start"), item.get("end")
            if start is None or end is None:
                raise ConfigError(f"{path}: '{block}' entry {item!r} needs both start and end")
            if int(start) < 1:
                raise ConfigError(f"{path}: '{block}' uses 1-based numbering; got start={start}")
            if int(end) < int(start):
                raise ConfigError(f"{path}: '{block}' entry {item!r} has end before start")


def load(
    config_path: str | os.PathLike[str],
    *,
    repo_root: str | os.PathLike[str] | None = None,
) -> Config:
    """Load and validate a target configuration."""
    path = Path(config_path).resolve()
    if not path.exists():
        raise ConfigError(f"no such config file: {path}")

    try:
        raw = yaml.safe_load(path.read_text()) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path}: not valid YAML -- {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: top level must be a mapping")

    name = (raw.get("project") or {}).get("name")
    if not name or not _SLUG.match(str(name)):
        raise ConfigError(
            f"{path}: project.name must be a lowercase slug (letters, digits, '_', '-'); "
            f"got {name!r}. It is used in output paths."
        )

    _check_weights(raw, path)
    _check_ranges(raw, path)

    root = Path(repo_root).resolve() if repo_root else Path(__file__).resolve().parents[2]
    return Config(raw=raw, path=path, _repo_root=root)
