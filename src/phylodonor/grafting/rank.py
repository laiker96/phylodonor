"""Composite scoring, filtering and output for the grafting stage.

Separated from :mod:`phylodonor.grafting.score` on purpose: scoring is
expensive and deterministic, ranking is cheap and full of judgement calls you
will want to revisit. Keeping them apart means re-weighting does not mean
re-running the language model.

What this module is optimising for
----------------------------------
This stage discards grafts that cannot work. It does not identify the best
one -- only expression and measurement do that. That makes **false negatives
the error that matters**: a good design silently filtered out is gone, while a
mediocre design that survives costs one well in a plate.

Three consequences, all of which show up in the code below:

- Filters are permissive by default, and a filter is only worth adding if its
  failure mode is well understood.
- Rejections are recorded with reasons, never dropped, and
  :func:`write_designs` emits a per-filter tally so an over-aggressive rule is
  visible rather than invisible.
- ``keep_top`` defaults to null. Truncating at a fixed rank assumes the
  ordering among survivors is trustworthy; none of the score terms is
  calibrated against measured stability, so it is not. The real cut is however
  many constructs you can express.

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
    """Score, filter and sort. Passing designs first, best composite first.

    Rejected designs stay in the returned list with ``passed_filters=False``
    so that callers can count and inspect them. Truncation to ``keep_top``, if
    set, happens here and applies only to passing designs.
    """
    raise NotImplementedError


def rejection_summary(ranked: Sequence[ScoredDesign]) -> dict[str, int]:
    """Count rejections per filter reason, plus totals.

    The diagnostic that matters for a tool whose job is discarding designs: if
    one rule is responsible for removing most of the candidate set, that is
    far more likely to be a bad threshold than a genuinely bad set of designs.
    Printed at the end of every run and written into the manifest.
    """
    raise NotImplementedError


def write_designs(ranked: Sequence[ScoredDesign], cfg: Config) -> Path:
    """Write ``designs.tsv``, ``designs.fa`` and ``manifest.json``.

    ``designs.tsv`` contains every design, passing and rejected, with all
    score terms and the rejection reasons -- it is the audit trail for what
    was thrown away. ``designs.fa`` contains only the survivors.

    The manifest records the config hash, the resolved language-model name,
    the tool versions and the :func:`rejection_summary`. A ranking without
    that provenance cannot be reproduced and should not be published.
    """
    raise NotImplementedError
