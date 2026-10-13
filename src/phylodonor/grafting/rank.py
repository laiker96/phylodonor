"""Composite scoring, filtering and output for the grafting stage.

Separated from :mod:`phylodonor.grafting.score` on purpose: scoring is
expensive and deterministic, ranking is cheap and full of judgement calls you
will want to revisit. Keeping them apart means re-weighting does not mean
re-running the language model.

TODO: implement. Signatures and contract are fixed.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from phylodonor.config import Config
from phylodonor.grafting.models import Design, LMScores, ScoredDesign

__all__ = ["composite_score", "apply_filters", "rank_designs", "write_designs"]


def composite_score(scores: LMScores, cfg: Config) -> float:
    """Weighted sum of the available score terms.

    Weights come from ``grafting.score.weights`` and are validated to sum to
    1.0 at config load. Terms that are ``None`` (because the corresponding
    model was not enabled) are dropped and the remaining weights
    renormalised, so that disabling inverse folding does not silently deflate
    every composite score.
    """
    raise NotImplementedError


def apply_filters(design: Design, scores: LMScores, cfg: Config) -> tuple[bool, tuple[str, ...]]:
    """Hard filters from ``grafting.score.filters``.

    Returns ``(passed, notes)``. Rejections are recorded rather than dropped:
    a filter that quietly removes most of the candidate set is a bug, and you
    can only see it if the rejections survive to the output table.
    """
    raise NotImplementedError


def rank_designs(
    designs: Sequence[Design],
    scores: Sequence[LMScores],
    cfg: Config,
) -> list[ScoredDesign]:
    """Score, filter and sort. Passing designs first, best composite first."""
    raise NotImplementedError


def write_designs(ranked: Sequence[ScoredDesign], cfg: Config) -> Path:
    """Write ``designs.tsv``, ``designs.fa`` and ``manifest.json``.

    The manifest records the config hash, the resolved language-model name and
    the tool versions. A ranking without that provenance cannot be reproduced
    and should not be published.
    """
    raise NotImplementedError
