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

## Two different questions

The toolkit computes both of the following, because they are not interchangeable and
treating them as such is the main way this stage goes wrong.

### Plausibility — likelihood

Score the pseudo-log-likelihood (PLL) of the chimeric sequence under a masked language
model and subtract the receptor's. Negative means the model finds the design less natural
than wild type.

This is the better-grounded signal: masked-token probability is exactly what the model was
trained to estimate, and it needs no reference set, no assumptions and no tuning. It is
very good at the thing you most need — rejecting designs that are quietly broken.

Two variants:

| `pll_method` | Cost | Use when |
|---|---|---|
| `masked` | one forward pass **per position** | default; required if you trust the junction term |
| `wt_marginal` | one forward pass **per design** | coarse first pass over a very large candidate set |

`wt_marginal` reads every position from a single unmasked pass. It correlates well with
`masked` for ranking overall, but it is biased precisely where local context changed — at
the junctions — which is where you most need accuracy. Use it to triage thousands down to
hundreds, then re-score with `masked`.

**The junction term is the one that earns its keep.** A graft can leave global PLL almost
untouched while producing two implausible seams, because the deltas average out over a
400-residue protein. Restricting the same per-position log probabilities to a window around
each seam surfaces exactly the failure mode that kills designs experimentally. If you
compute one thing, compute this.

### Direction — embedding distance

Embed the chimera, embed a reference set (your thermophilic donors, or the tips of the
thermophilic clades from the `thermal` stage), and measure the distance between them.

This asks a question likelihood cannot: *has this design actually moved toward the
thermophilic regime?* Likelihood is objective-blind — a design can be perfectly plausible
and no more thermostable than wild type. Embedding distance encodes the goal directly, and
costs one forward pass per design.

Its weakness is that it is only as meaningful as the reference set, and embedding proximity
has no calibrated relationship to stability *or* foldability. A design can drift toward the
thermophile centroid while being unfoldable, and nothing in the metric notices.

### How to combine them

Use likelihood to **reject** and embedding distance to **order**:

1. Hard-filter on `max_pll_drop` — remove designs the model finds implausible.
2. Hard-filter on the junction term — remove designs with bad seams.
3. Rank survivors by the composite, with embedding distance carrying real weight.

Do not ask a single weighted sum to do both jobs. A design with a catastrophic junction and
an excellent embedding distance will score respectably in a naive composite, and it will
fail in the lab.

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
