# External tools and licensing

This repository is MIT licensed. The tools it *invokes* are not, and some of them restrict
what you may do with the results. None are bundled; all are installed by you.

Read this before using the pipeline in a commercial setting.

## Licence-restricted — check before commercial use

| Tool | Used by | Restriction |
|---|---|---|
| **Amber / pmemd** | `dynamics` (Amber engine) | Commercial licence required. The source archive must be obtained from the Amber distributors and supplied as a local path; it is never tracked here. The OpenMM engine is a free alternative covering the same protocols. |
| **PyRosetta** | optional scoring | Free for non-commercial and academic use only. A commercial licence is separate. |
| **CamSol** | `scoring.solubility` | The published method's fitted coefficients are proprietary and are **not** included. Any local reimplementation is therefore *not* equivalent to the published predictor and must not be described as CamSol results. Use the web server for publishable numbers. |
| **vmd-python** | optional trajectory analysis | VMD's licence is free for non-commercial use but is not an OSI open-source licence. |

## Free, but cite them

These are the load-bearing dependencies. If you publish results, cite them — the pipeline
is a wrapper around their work.

| Tool | Used by | Licence |
|---|---|---|
| MMseqs2 | `homologs` | MIT |
| HMMER | `phylogeny` | BSD-3 |
| trimAl | `phylogeny` | GPL-3 |
| IQ-TREE | `phylogeny` | GPL-2 |
| OpenMM | `dynamics` | MIT / LGPL |
| ESM models | `grafting`, `scoring` | MIT (weights: check each release) |
| US-align / TM-align | `donors` | free for academic use; author permission for commercial |
| Foldseek | `donors` | GPL-3 |
| Arpeggio / GetContacts | `holo` | GPL-3 / MIT |
| DSSP | `grafting` | BSD-2 |
| Aggrescan3D | `scoring.aggregation` | free web server + academic licence |

GPL tools invoked as **separate processes** do not make this toolkit a derivative work. Do
not link them into the package.

## Model weights

Language-model weights carry their own terms, separate from the code that loads them. Check
the licence of the specific checkpoint you configure in
`grafting.language_model.sequence_model` — within a single model family, different releases
have shipped under different terms, including some that are non-commercial.

## Reference data

Growth-temperature datasets and sequence databases are third-party data with their own
citation and redistribution terms. None are vendored here; see
[data_sources.md](data_sources.md) for how to obtain them and what each one requires.

## Structures and sequences you supply

Anything you point the config at — your target sequence, predicted structures, curated
design tables — stays yours and stays out of this repository. `data/targets/`,
`data/reference/` and `structures/` are gitignored except for their `.gitkeep` files,
deliberately: it should be difficult to commit your own inputs by accident.
