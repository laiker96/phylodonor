"""Data structures shared across the grafting stage.

These are the stage's vocabulary. Keeping them in one place means the
enumeration, scoring and ranking steps can be developed and tested
independently, and that a design can be serialised to TSV and read back
without losing provenance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

__all__ = ["Segment", "GraftSpec", "Design", "LMScores", "ScoredDesign"]


@dataclass(frozen=True, slots=True)
class Segment:
    """A contiguous residue range in one protein, 1-based and inclusive."""

    source: str
    """Accession or identifier of the protein this range belongs to."""
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 1:
            raise ValueError(f"{self}: coordinates are 1-based; got start={self.start}")
        if self.end < self.start:
            raise ValueError(f"{self}: end precedes start")

    def __len__(self) -> int:
        return self.end - self.start + 1

    def label(self) -> str:
        return f"{self.source}:{self.start}-{self.end}"


@dataclass(frozen=True, slots=True)
class GraftSpec:
    """One substitution: a receptor segment replaced by a donor segment.

    The two segments need not be the same length -- that is the point of
    grafting a different loop or helix -- but the length difference is
    constrained by ``grafting.candidates.max_length_delta``.
    """

    receptor: Segment
    donor: Segment
    rationale: str = ""
    """Short human-readable reason this pair was proposed. Carried through to
    the output table so a reviewer can tell a structural graft from a
    disulfide design without re-deriving it."""

    @property
    def length_delta(self) -> int:
        return len(self.donor) - len(self.receptor)

    def label(self) -> str:
        return f"{self.receptor.start}-{self.receptor.end}<-{self.donor.label()}"


@dataclass(frozen=True, slots=True)
class Design:
    """A candidate chimera: the receptor with one or more grafts applied."""

    design_id: str
    grafts: tuple[GraftSpec, ...]
    sequence: str
    """Full-length chimeric sequence, grafts already applied."""
    junctions: tuple[int, ...] = ()
    """Positions in ``sequence`` (1-based) of each seam. Two per graft, unless
    adjacent grafts merge. Scoring windows are centred on these."""

    def __post_init__(self) -> None:
        if not self.grafts:
            raise ValueError(f"{self.design_id}: a design needs at least one graft")
        if not self.sequence:
            raise ValueError(f"{self.design_id}: empty sequence")


@dataclass(frozen=True, slots=True)
class LMScores:
    """Language-model terms for one design, all as deltas against the receptor.

    Zero means "indistinguishable from wild type". Negative means the model
    considers the design less plausible. Values are per-residue averages so
    that designs of different lengths stay comparable.
    """

    embedding_distance: float | None = None
    """Movement of the design's mean embedding toward a thermophilic reference
    set, relative to the receptor. Positive means the design sits closer to
    the reference than wild type does.

    This is the default ranking signal, because it is the one that has been
    run end to end here with structure-prediction follow-up. It encodes the
    objective directly -- likelihood does not know you want thermostability --
    at one forward pass per design. Its limitation is that it is only as
    meaningful as the reference set, and embedding proximity has no calibrated
    relationship to stability or foldability."""
    donor_conservation: float | None = None
    """How typical the donor segment is within its own thermophilic clade. A
    graft taken from an idiosyncratic donor loop is a worse bet than one taken
    from a region conserved across thermophiles."""

    delta_pll: float | None = None
    """Whole-sequence pseudo-log-likelihood, chimera minus receptor.

    Better motivated in principle than the embedding term -- masked-token
    probability is what the model was trained to estimate -- but unvalidated
    in this pipeline, so it carries no weight by default. See
    ``docs/grafting.md``."""
    junction_pll: float | None = None
    """Mean delta PLL within the junction windows only. Expected to be the
    most discriminating likelihood term, since a graft can be globally
    unremarkable while having two implausible seams, and seams are where
    designs fail. Also unvalidated here."""
    inverse_folding: float | None = None
    """Sequence-given-backbone log-likelihood delta, if an inverse-folding
    model was enabled. Requires a predicted structure per design."""

    model_name: str = ""
    pll_method: Literal["masked", "wt_marginal", ""] = ""


@dataclass(frozen=True, slots=True)
class ScoredDesign:
    """A design with its language-model terms and composite rank."""

    design: Design
    scores: LMScores
    composite: float
    passed_filters: bool
    filter_notes: tuple[str, ...] = field(default_factory=tuple)
    """Why a design was rejected, if it was. Recorded rather than discarded --
    a filter that silently eats most of the candidate set is a bug you want to
    be able to see."""
