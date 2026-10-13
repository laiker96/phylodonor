# Roadmap

The architecture, configuration schema and stage contracts are settled. Stages are being
implemented in dependency order, easiest first, so that a runnable slice exists early.

## Order

| # | Stage | Status | Notes |
|---|---|---|---|
| 1 | `holo` | in progress | Cofactor transfer from donor crystals, protonation, restrained minimisation. Few dependencies, so it lands first. |
| 2 | folding helper | in progress | Sequence → structure, used by `holo` and `grafting`. A thin wrapper over an existing predictor rather than its own stage. |
| 3 | `dynamics` | planned | Constant-temperature and heating-ramp protocols. Single implementation behind `dynamics.engine`, so that OpenMM and Amber share one analysis path instead of two. |
| 4 | `donors` | planned | Structural superposition and per-domain donor ranking. Needs the region-map abstraction below. |
| 5 | `grafting` | planned | Candidate enumeration and language-model triage. Interfaces are already fixed in `src/phylodonor/grafting/`; see [grafting.md](grafting.md). |
| 6 | `homologs`, `phylogeny`, `thermal` | planned | The reconnaissance half. Largest piece of work, and the one that most needs the region map. |
| 7 | `scoring`, `viz` | planned | Predictors and figures. Useful throughout but blocking nothing. |

Stages 1–4 form a coherent tool on their own: fold → build a cofactor-bound model → find
thermophilic donors → test them in simulation. That is the first milestone.

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
