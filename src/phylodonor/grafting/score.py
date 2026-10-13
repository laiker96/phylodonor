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

**Likelihood** (``delta_pll``, ``junction_pll``). Ask the model how probable
the chimeric sequence is, and compare against the unmodified receptor. This
measures *plausibility*: would evolution plausibly have produced this? It is
the better-grounded choice, since masked-token likelihood is exactly what the
model was trained to estimate, and it needs no reference set. Its weakness is
cost -- the masked-marginal estimator is one forward pass per position.

**Embedding similarity** (``embedding_distance``). Embed the chimera and
measure distance to a reference set, typically the thermophilic donors. This
asks *does this look more like a thermophile than wild type did?* It is far
cheaper (one forward pass per design) and it encodes the thermostability
objective directly, which likelihood does not. Its weakness is that it is only
as meaningful as the reference set, and embedding proximity has no calibrated
relationship to either stability or foldability -- a design can drift toward
the thermophile centroid while being unfoldable.

Neither is a substitute for simulation or experiment. Use likelihood to
*reject* implausible designs and embedding distance to *order* the survivors
by how far they move toward the thermophilic regime, rather than asking one
number to do both jobs.

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
        Donor or thermophile-clade sequences, used only for the optional
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
