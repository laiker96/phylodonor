# Ranking grafts with a protein language model

This is the design stage, and the one where the toolkit makes the most consequential
choices. Everything upstream is reconnaissance; this is where candidate proteins get
proposed and triaged.

## The problem

You have a receptor (your target) and a handful of thermophilic donors with per-domain
structural superpositions. For each donor there are many structurally plausible segments
you could transplant, each with some freedom in exactly where the junctions fall, and
designs may combine more than one graft. A modest enumeration — four donors, six candidate
regions each, a few junction offsets, pairs allowed — reaches the low thousands.

You can afford to fold tens and simulate a handful. So something has to cut three orders of
magnitude, cheaply, without discarding the good designs.

Geometric filters alone are not enough. They tell you a segment fits; they do not tell you
whether the resulting sequence is something that folds and expresses. That is what a protein
language model is for.

## Validation status — read this before changing the weights

The two scoring approaches below are not on equal footing in this pipeline, and the
defaults encode that rather than encoding which one is prettier in theory.

| Path | Status |
|---|---|
| **Embedding distance** | Run end to end on a real target. Ranked designs were folded with a structure predictor and the results held up. This is the default. |
| **Likelihood (PLL)** | Implemented, **not run**. Carries zero weight by default, and `filters.max_pll_drop` defaults to `null`. |

Likelihood is the better-motivated choice on paper, for the reasons below. It is still not
the default, because a hard filter on an unvalidated signal silently discards designs you
have no evidence are bad — and silent discarding is the failure mode you cannot detect from
the output. Switching is a config change; do it against designs whose structural or
experimental outcome you already know.

## Two different questions

The toolkit computes both, because they are not interchangeable and treating them as such
is the main way this stage goes wrong.

### Direction — embedding distance (default)

Embed the chimera, embed a reference set, and measure how far the design moves toward it
relative to wild type. `embedding_reference` selects the set: `donors` uses the ranked
donors, `thermophile_clades` uses every tip the `thermal` stage called thermophilic — a
larger and less arbitrary reference, and the better choice when the tree supports it.

This asks a question likelihood cannot: *has this design actually moved toward the
thermophilic regime?* Likelihood is objective-blind — a design can be perfectly plausible
and no more thermostable than wild type. Embedding distance encodes the goal directly, at
one forward pass per design.

Its weakness is that it is only as meaningful as the reference set, and embedding proximity
has no calibrated relationship to stability *or* foldability. A design can drift toward the
thermophile centroid while being unfoldable, and nothing in the metric notices. **This is
why the structure-prediction gate below is not optional.**

### Plausibility — likelihood

Score the pseudo-log-likelihood of the chimeric sequence under a masked language model and
subtract the receptor's. Negative means the model finds the design less natural than wild
type.

Better grounded in principle: masked-token probability is exactly what the model was
trained to estimate, it needs no reference set, and it should be good at the thing the
embedding term cannot do — rejecting designs that are quietly broken.

Two variants:

| `pll_method` | Cost | Use when |
|---|---|---|
| `masked` | one forward pass **per position** | required if you trust the junction term |
| `wt_marginal` | one forward pass **per design** | coarse first pass over a very large candidate set |

`wt_marginal` reads every position from a single unmasked pass. It correlates well with
`masked` for overall ranking, but it is biased precisely where local context changed — at
the junctions — which is where you most need accuracy.

**If you validate one likelihood term, validate the junction term.** A graft can leave
global PLL almost untouched while producing two implausible seams, because the deltas
average out over a 400-residue protein. Restricting the same per-position log probabilities
to a window around each seam targets exactly the failure mode that kills designs. That is
the argument for it; it is not yet a result.

### The combination to aim at

Likelihood to **reject**, embedding distance to **order**:

1. Hard-filter on `max_pll_drop` — remove designs the model finds implausible.
2. Hard-filter on the junction term — remove designs with bad seams.
3. Rank survivors by embedding distance and donor conservation.

Do not ask a single weighted sum to do both jobs. A design with a catastrophic junction and
an excellent embedding distance scores respectably in a naive composite and fails in the
lab. Note this is the *target* design, not the validated one — step 3 is what currently
runs, with steps 1 and 2 available and unproven.

## The structure-prediction gate

Because the default ranking signal says nothing about foldability, the survivors get folded
and checked. `validate.junction_confidence` gates on predicted confidence **at the seams**
rather than over the whole model: a broken graft typically yields a confident global fold
with one disordered junction, and a mean pLDDT across 400 residues averages that away.

What this establishes, precisely: it rules out structural incoherence. It is not evidence
of stability, and structure predictors are permissive about chimeras — they will
confidently fold sequences that never express. Passing the gate means a design is worth
spending simulation time on, nothing more.

## Optional: inverse folding

An inverse-folding model (ESM-IF, ProteinMPNN) scores the sequence *conditioned on a
backbone*, which sharpens the question from "is this a plausible protein" to "is this a
plausible protein for this fold". That is closer to what you actually want to know.

The cost is a structure per design, so it is off by default. Enable it once the candidate
set is down to tens — at which point it is the most informative term available short of
simulation.

## What this stage does not tell you

Likelihood is not stability. The model has no representation of your assay temperature,
your expression host, your cofactor, or your buffer. Published correlations between
language-model scores and measured ΔΔG are real but modest, and they are weakest for the
multi-residue structural changes this pipeline generates — most benchmarks are built on
single point mutations.

So: the ordering *within* your top ten is not trustworthy. Treat this stage as a filter
that gets you from thousands of designs to a shortlist worth folding, and let dynamics and
experiment decide among the shortlist. Designs that fail the filters are almost always
genuinely bad; designs that pass it are merely *not obviously bad*, which is all you need
from a triage step.

## Configuration

See the `grafting:` block in [`configs/template.yaml`](../configs/template.yaml). The
fields that matter most:

- `candidates.exclude_positions` — catalytic residues, cofactor contacts, interfaces. The
  single most important field in the whole config. Unconstrained enumeration proposes
  grafts through the active site, and the language model will not object: an active site
  replaced by plausible-looking surface loop scores *better* than wild type, because active
  sites are strained and unusual. The model cannot protect you here; only this list can.
- `candidates.junction_rules` — keeps seams in exposed loops where the backbones already
  superpose.
- `score.weights` — validated to sum to 1.0 at config load. Terms that are unavailable are
  dropped and the rest renormalised.
- `score.filters.forbid_new_cysteines` — set `true` unless you are deliberately designing
  disulfides. A free cysteine introduced by accident is a recurring source of aggregation
  and heterogeneous expression.

## References

See [REFERENCES.md](../REFERENCES.md) for the masked-marginal estimator, the ESM model
family, inverse folding, and benchmarks on language-model stability prediction.
