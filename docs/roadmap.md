# Roadmap

The architecture, configuration schema and stage contracts are settled. Stages are being
implemented in dependency order, easiest first, so that a runnable slice exists early.

## What is required, and what is not

The pipeline this was extracted from had nine stages. Five are required here. The cuts were
made by tracing what actually consumes each output, not by taste.

| Dropped from the critical path | Why |
|---|---|
| Phylogenetic tree | Contributes nothing to donor ranking; the donor score is per-domain structural similarity, coverage and TM-score, with no phylogenetic term. Tree distances were computed and then used only as a heatmap sort key, in a module nothing imported. |
| Molecular dynamics | No stage reads it. See below. |
| Solubility / aggregation scoring | Fed only a terminal analysis step. Gates no design. |
| Sequence logos, conservation profiles, interactive trees, clade differential analysis | Figures. Nothing reads them. |

In the source pipeline, six of twenty-six steps in the largest sub-pipeline were terminal
leaves whose output nothing consumed, and the whole sub-pipeline exported exactly one
artifact to the rest of the project.

### The one thing the tree really did

It cleaned growth-temperature labels. Clade calling required several annotated tips that all
agreed on regime within a supported subtree, which discards singleton and suspect
annotations — and since growth-temperature datasets join by taxon id, those are precisely
the errors worth discarding.

That function is kept, by a cheaper mechanism: `thermal` calls regimes over the **sequence
clusters** that `homologs` already computes for representative selection. Same logic —
agreement among related sequences — without the alignment, trimming and tree inference.

What is deliberately given up: imputing a regime across deeper phylogenetic distance. The
original could promote an unannotated sequence into the thermophile pool on clade
membership; `impute_unlabelled` defaults to false here, because that is a real claim to make
about a sequence nothing is known about.

### A note on where the design insight came from

Worth recording, because it is the real architectural difference. In the source pipeline the
graft targets were hand-curated constants, and the analysis steps that informed them passed
their conclusions through a *person*, not a file — no code path connects them. Those steps
are data leaves, but they were not useless; they were how someone decided where to graft.

Here that step is mechanical: enumerate every structurally plausible graft under explicit
constraints, then triage with a language model and a structure-prediction gate. The analysis
steps can go not because the thinking was worthless but because candidates are now generated
and filtered rather than chosen.

## Simulation is not on the critical path

Worth stating first, because it shapes everything below.

**Minimisation and molecular dynamics are different things.** Minimisation is a
deterministic descent to a nearby energy minimum: seconds to minutes, no trajectory. MD
integrates motion over time: GPU-days for a panel of designs. The original codebase treated
them as one concern; they are not.

Minimisation is in `holo`, and even there it is optional — the contact table that `grafting`
consumes does not need it, since relaxation shifts atoms by tenths of an angstrom and rarely
changes which residues are in contact.

**MD is opt-in and off the critical path.** What it was providing:

| Was used for | Replaced by |
|---|---|
| A relaxed holo model | Restrained minimisation, which is not MD |
| Per-residue RMSF, for junction placement | DSSP secondary structure, relative SASA, and per-residue predicted confidence (pLDDT), which is itself a disorder signal |
| Heating-ramp apparent unfolding, for ranking | Nothing — see below |

That last row is the substantive decision. Short MD is an unreliable estimator of ΔTm: 100 ns
of a 400-residue protein samples essentially none of the unfolding landscape, and a ramp
protocol's apparent transition depends heavily on heating rate with no clean mapping to a
measured melting temperature. Since this pipeline's job is to *discard* grafts that cannot
work — and a bad graft is already caught by the language model, the clash detector and the
structure-prediction gate, all orders of magnitude cheaper — MD would only offer to rank,
which is the thing not worth trusting. Experiment decides that.

Consequences: the Amber dependency disappears from the critical path, so the only
licence-encumbered hard dependency is gone; `holo` is standalone; and the whole pipeline runs
on a laptop plus a folding server.

MD stays in the tree as an opt-in stage. It earns its place for a specific mechanistic
question ("does this loop open before the core melts?"), or retrospectively if measured Tm
values ever arrive and the correlation is worth checking.

## Order

| # | Stage | Status | Notes |
|---|---|---|---|
| 1 | `holo` | **done** | Cofactor transfer, optional restrained minimisation, contact annotation. |
| 2 | folding helper | next | Sequence → structure, used by `holo` and `grafting`. A thin wrapper over an existing predictor rather than its own stage. |
| 3 | `donors` | planned | Structural superposition and per-domain donor ranking. Needs the region-map abstraction below. |
| 4 | `grafting` | planned | Candidate enumeration and language-model triage. Interfaces are fixed in `src/phylodonor/grafting/`; see [grafting.md](grafting.md). |
| 5 | `homologs`, `thermal` | planned | The reconnaissance half: search, subfamily filter, cluster, align, then regimes over clusters. Substantially smaller than the source pipeline now that the tree is gone. |
| — | `scoring`, `viz`, `phylogeny`, `dynamics` | opt-in | Nothing depends on any of them. Port if a question needs one. |

Stages 1–4 form a coherent tool on their own: fold → build a cofactor-bound model → find
thermophilic donors → triage grafts. That is the first milestone, and it needs no GPU.

## The region-map abstraction

The one piece of design debt worth naming up front. Several stages need to know a protein's
domain decomposition: donor scoring is per domain, sequence logos are per domain, alignment
refinement is per inter-domain linker, and output tables are keyed by domain name.

The wrong way to do this — and the easy way — is a module of constants describing one
protein's domains, imported wherever needed. It works until the second protein, at which
point the domain names have leaked into column headers, scoring weights and filenames, and
every one of those is a breaking change.

So: domains come from `domains.definitions` in the config, stages receive them as data, and
anything derived from them (column names, weights, output keys) is generated rather than
written down. `donors` and the MSA stages are both blocked on this, which is why they are
late in the order despite being scientifically central.

## Design rules

These apply to every stage and are worth stating once:

- **No protein-specific value in `src/`.** Not an accession, a residue number, an organism
  name, a sequence or an input path. See [CONTRIBUTING.md](../CONTRIBUTING.md).
- **Coordinates are 1-based and inclusive**, converted once at any library boundary.
- **Receptor numbering is canonical**; donor coordinates are always derived through an
  explicit structural mapping.
- **Every stage emits a manifest** with the config hash and the versions of every external
  tool and model used.
- **Licence-encumbered tools stay optional.** Nothing that restricts commercial use may
  become a core dependency. See [external_tools.md](external_tools.md).

## Known gaps

- Two of the planned MSA steps consume structural annotation tables (cofactor contacts,
  per-pocket conservation) that no stage currently produces. Either `holo` generates them or
  those steps become optional; decide before stage 6.
- The `dynamics` analysis code should be written once against a trajectory abstraction. Two
  engine-specific copies of the same RMSF and timeseries analysis is a maintenance trap.
- No end-to-end integration test yet. The demo config exists for this; it needs a tiny
  fixture database so CI can run the whole pipeline on a handful of sequences.
