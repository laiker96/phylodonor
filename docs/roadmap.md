# Roadmap

The architecture, configuration schema and stage contracts are settled. Stages are being
implemented in dependency order, easiest first, so that a runnable slice exists early.

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
| 5 | `homologs`, `phylogeny`, `thermal` | planned | The reconnaissance half. Largest piece of work, and the one that most needs the region map. |
| 6 | `scoring`, `viz` | planned | Predictors and figures. Useful throughout but blocking nothing. |
| — | `dynamics` | opt-in | Not required by any other stage. Port when a question needs it. |

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
