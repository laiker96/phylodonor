# phylodonor

Phylogeny-guided donor grafting for protein thermostability, triaged by a protein language model.

> **Status: scaffold.** The architecture, configuration schema and stage contracts are in
> place; stage implementations are being ported in. See [docs/roadmap.md](docs/roadmap.md).

## What it does

Given an enzyme you want to make more thermostable, and nothing else:

```
homologs   search a sequence database, filter to the right subfamily, cluster,
           and align the representatives
thermal    annotate each homolog's source organism with growth temperature and
           call psychro- / meso- / thermophilic regimes over the clusters
donors     rank thermophilic homologs by how well their domains superpose on
           your target -- these are your donors of stabilising structure
holo       build cofactor-bound models of target and donors; annotate contacts
grafting   enumerate chimeric designs that transplant donor elements into the
           target, and triage them with a protein language model
                                 ↓
           a shortlist to express

optional, nothing above depends on them:
scoring    conservation, solubility and aggregation predictors
viz        logos and plots
phylogeny  alignment trimming and tree inference
dynamics   molecular dynamics
```

**Five required stages, and no GPU beyond the language model.** Neither a phylogenetic
tree nor molecular dynamics is required, and both were in the pipeline this was extracted
from:

- A **tree** contributes nothing to donor ranking — the donor score has no phylogenetic
  term. Its one real job was discarding noisy growth-temperature labels, which `thermal`
  now does by requiring agreement within a sequence cluster: the same mechanism, reusing
  clustering the pipeline already performs, without an alignment-trim-tree chain.
- **MD** was providing a relaxed model (that is minimisation, not MD), per-residue
  flexibility for junction placement (DSSP, SASA and predicted confidence cover it), and a
  heating-ramp stability ranking — which short MD estimates unreliably, and which this
  pipeline deliberately does not try to do.

External tools reduce to mmseqs2, a structural aligner, DSSP and ESM. Nothing
licence-encumbered. [docs/roadmap.md](docs/roadmap.md) has the reasoning.

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
