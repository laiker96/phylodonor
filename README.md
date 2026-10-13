# phylodonor

Phylogeny-guided donor grafting for protein thermostability, triaged by a protein language model.

> **Status: scaffold.** The architecture, configuration schema and stage contracts are in
> place; stage implementations are being ported in. See [docs/roadmap.md](docs/roadmap.md).

## What it does

Given an enzyme you want to make more thermostable, and nothing else:

```
homologs   find family members in a sequence database, filter to the right subfamily
phylogeny  align, trim, infer a tree
thermal    annotate each homolog's source organism with growth temperature,
           call clades that are consistently psychro- / meso- / thermophilic
donors     rank thermophilic homologs by how well their domains superpose on
           your target -- these are your donors of stabilising structure
holo       build cofactor-bound models of target and donors; annotate contacts
grafting   enumerate chimeric designs that transplant donor elements into the
           target, and rank them with a protein language model
scoring    conservation, solubility and aggregation predictors
dynamics   optional: simulate survivors, for a specific mechanistic question
```

No GPU and no molecular dynamics are required. Minimisation in `holo` is a seconds-long
restrained relaxation and is itself optional; MD is an opt-in stage that nothing else
depends on. [docs/roadmap.md](docs/roadmap.md) explains why — briefly, short MD is an
unreliable estimator of melting temperature, and this pipeline's job is to discard grafts
that cannot work rather than to rank the ones that might.

The premise: a thermophilic relative of your enzyme has already solved the stability
problem, and the structural elements responsible can be transplanted. The tree tells you
which relatives are genuinely thermophilic rather than incidentally annotated; the
structural comparison tells you which of their elements correspond to yours; the language
model tells you which of the many possible transplants are worth building.

### Why a language model in the loop

Enumerating grafts across a handful of donors produces hundreds to thousands of candidate
chimeras. Folding and simulating all of them is not affordable, and geometric filters alone
do not distinguish a graft that reads like a real protein from one that merely fits.

A protein language model gives a cheap sequence-only signal. Two are implemented, and they
answer different questions:

- **Direction** — how far the design's embedding moves toward a thermophilic reference set,
  relative to wild type. This encodes the actual objective, which likelihood does not.
- **Plausibility** — pseudo-log-likelihood of the chimera versus the receptor, globally and
  in windows around each junction. A graft can leave global likelihood almost untouched
  while producing two implausible seams, and seams are where designs fail.

Direction is the default, because it is the path that has been run end to end with
structure-prediction follow-up. The likelihood terms are better motivated in principle but
unvalidated here, so they carry no weight by default — switching is a config change, not a
code change. [docs/grafting.md](docs/grafting.md) lays out the status of each.

Neither signal is calibrated against stability. This stage is a filter that takes thousands
of designs down to tens worth simulating, not an oracle that picks the winner.

## Install

```bash
git clone https://github.com/laiker96/phylodonor.git
cd phylodonor
pip install -e ".[dev]"            # core + tooling
pip install -e ".[lm,structure]"   # add language models and structure handling
```

Heavy and licence-encumbered dependencies are deliberately optional extras: `lm` pulls
PyTorch, `dynamics` pulls OpenMM. Some tools the pipeline can call (Amber/pmemd, PyRosetta,
CamSol) are separately licensed and are **not** bundled — see
[docs/external_tools.md](docs/external_tools.md).

## Use

Every stage reads one YAML config. Nothing in `src/` hardcodes a protein.

```bash
cp configs/template.yaml configs/my_target.yaml
$EDITOR configs/my_target.yaml

phylodonor homologs search --config configs/my_target.yaml --dry-run   # validate first
phylodonor homologs search --config configs/my_target.yaml
phylodonor phylogeny tree  --config configs/my_target.yaml
phylodonor grafting  rank  --config configs/my_target.yaml
```

`configs/template.yaml` documents every field. The two that most affect results:

- `domains.definitions` — the domain decomposition of your target. Donors are scored per
  domain, and a wrong decomposition produces confidently wrong donor rankings.
- `grafting.candidates.exclude_positions` — catalytic residues, cofactor contacts, binding
  interfaces. Unconstrained enumeration will cheerfully propose grafts through the active
  site.

Reference databases (UniRef, growth-temperature datasets) are **not** downloaded for you;
point the config at local copies. See [docs/data_sources.md](docs/data_sources.md).

## Keeping your target data out of the repository

Protein engineering projects are usually confidential well before they are published, and
this pipeline is designed so that nothing target-specific ever needs to live in the code:
sequences, structures, accessions, domain boundaries and design tables all come from your
config and your own input directories, which are gitignored by default.

`tools/check_denylist.py` makes that a guarantee rather than an intention. It runs as a
pre-commit hook and in CI, over both the working tree and the full git history, and fails
on personal absolute paths, credential shapes, REPL history files and any pattern you give
it. Your own identifiers go in a gitignored `.private-patterns` file — a pattern list naming
what it protects does not belong in a public tree — or in the `PRIVATE_PATTERNS` CI secret:

```bash
pre-commit install                      # do this before your first commit
python tools/check_denylist.py --all
```

A pattern caught at commit time is a non-event; the same pattern caught after a push means
rewriting history.

## Citing

See [CITATION.cff](CITATION.cff). Methods this pipeline builds on are credited in
[REFERENCES.md](REFERENCES.md); if you publish results obtained with it, cite the
underlying tools (MMseqs2, HMMER, IQ-TREE, ESM, OpenMM) as well.

## Licence

MIT, for the code in this repository. Externally invoked tools carry their own licences,
some of which restrict commercial use.
