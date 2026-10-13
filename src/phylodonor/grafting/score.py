"""Language-model scoring of candidate grafts.

Scope
-----
This module turns a :class:`~phylodonor.grafting.models.Design` into an
:class:`~phylodonor.grafting.models.LMScores`. It does not decide anything --
filtering and ranking live in :mod:`phylodonor.grafting.rank`.

Scoring approaches
------------------
There are two defensible ways to score a chimera with a protein language
model, and they answer different questions.

**Embedding similarity** (``embedding_distance``) -- the default. Embed the
chimera and measure how far it moves toward a thermophilic reference set
relative to the receptor. This asks *does this look more like a thermophile
than wild type did?* It is cheap (one forward pass per design) and it encodes
the thermostability objective directly, which likelihood does not. Its
weakness is that it is only as meaningful as the reference set, and embedding
proximity has no calibrated relationship to stability or foldability -- a
design can drift toward the thermophile centroid while being unfoldable.

**Likelihood** (``delta_pll``, ``junction_pll``). Ask the model how probable
the chimeric sequence is and compare against the receptor. This measures
*plausibility* rather than direction, and it is better motivated in
principle: masked-token probability is exactly what the model was trained to
estimate, it needs no reference set, and the junction term targets the
specific failure mode that kills grafts. Its costs are compute (the
masked-marginal estimator is one forward pass per position) and, here,
evidence -- see below.

Validation status
-----------------
Be aware of the asymmetry, because the defaults encode it. The embedding path
has been run end to end on a real target, with the ranked designs checked by
structure prediction, and it held up. The likelihood path has not been run.

So likelihood carries no weight in the default config despite being the
better-motivated choice, and ``filters.max_pll_drop`` defaults to ``null``: a
hard filter on an unvalidated signal silently discards designs there is no
evidence are bad. Both terms are implemented and switching the weights is a
config change, but validate against designs with known structural or
experimental outcomes before trusting them to rank.

Neither signal is a substitute for simulation or experiment. The eventual
sensible combination is likelihood to *reject* and embedding distance to
*order*, rather than one weighted sum doing both jobs -- but that is a
hypothesis about this pipeline, not a result from it.

Implementation notes
--------------------
``pll_method="masked"`` is the masked-marginal estimator of Meier et al.
(2021): for each position, mask it, run the model, and take the log
probability of the observed residue. ``pll_method="wt_marginal"`` reads all
positions from a single unmasked pass -- roughly L times cheaper, well
correlated for ranking, and biased wherever a graft changes local context,
which is precisely at the junctions. Prefer ``masked`` when the junction term
matters, which is most of the time.

All returned values are per-residue averages, so designs of differing length
remain comparable.

TODO: implement. The signatures and the contract below are fixed; the bodies
are the work.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Protocol

from phylodonor.config import Config
from phylodonor.grafting.models import Design, LMScores

__all__ = ["SequenceScorer", "score_designs", "load_scorer"]


class SequenceScorer(Protocol):
    """Minimal interface a language model must satisfy to be used here.

    Defined as a Protocol so that ESM-2, ESM-C, a future model, or a cheap
    stub for tests are all interchangeable without touching the ranking code.
    """

    name: str

    def log_probs(self, sequences: Sequence[str], *, masked: bool) -> list[list[float]]:
        """Per-position log probability of each observed residue.

        Returns one list of length ``len(seq)`` per input sequence.
        """
        ...

    def embed(self, sequences: Sequence[str]) -> list[list[float]]:
        """Mean-pooled residue embedding per sequence."""
        ...


def load_scorer(cfg: Config) -> SequenceScorer:
    """Instantiate the scorer named by ``grafting.language_model``.

    Resolves device (``auto`` -> cuda if available), loads weights, and
    records the resolved model name so it can be written into the stage
    manifest. Model identity belongs in provenance: a ranking is not
    reproducible without it.
    """
    raise NotImplementedError


def score_designs(
    designs: Iterable[Design],
    receptor_sequence: str,
    cfg: Config,
    *,
    scorer: SequenceScorer | None = None,
    reference_sequences: Sequence[str] | None = None,
) -> list[LMScores]:
    """Score each design against the receptor.

    Parameters
    ----------
    designs
        Candidates from :mod:`phylodonor.grafting.candidates`.
    receptor_sequence
        The unmodified target sequence. Every term is a delta against this, so
        it is scored once and reused.
    cfg
        Supplies ``grafting.language_model.*`` and ``grafting.score.junction_window``.
    scorer
        Injected for testing; otherwise built via :func:`load_scorer`.
    reference_sequences
        Donor or thermophilic-homolog sequences, used only for the optional
        ``embedding_distance`` term. Omit to skip it.

    Returns
    -------
    One :class:`LMScores` per design, in input order.

    Notes
    -----
    Order of operations matters for cost: score the receptor once, batch
    designs according to ``batch_size``, and compute the junction term from
    the same forward passes as the global term rather than re-running the
    model over windows.
    """
    raise NotImplementedError


def _junction_windows(design: Design, window: int, seq_len: int) -> list[tuple[int, int]]:
    """Half-open index ranges around each junction, clipped to the sequence.

    Overlapping windows from adjacent grafts are merged so that residues in
    the overlap are not counted twice in the mean.
    """
    raise NotImplementedError
