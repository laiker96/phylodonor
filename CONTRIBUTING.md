# Contributing

## Setup

```bash
pip install -e ".[dev]"
pre-commit install          # required -- see below
pytest
```

## The one rule

**No protein-specific value may appear in `src/`.** Not an accession, a residue number, an
organism name, a melting temperature, a sequence, or an input path. If a stage needs one, it
reads it from the config via `phylodonor.config`.

This is what makes the toolkit reusable across targets, and it is also what makes it safe to
work on a confidential project in a public repository — the guarantee that no target data is
in the package is only as good as this rule.

One pre-commit hook enforces the part that cannot be undone: it refuses any run of 25+
amino-acid letters, and structure or sequence files by extension. Install it before your
first commit — a push to a public remote is permanent.

```bash
pre-commit install
```

## Adding a stage or step

1. Add the step name to `STAGES` in `src/phylodonor/cli.py`.
2. Add required config keys to `_STAGE_REQUIREMENTS` in `src/phylodonor/config.py`, so a
   missing database path fails in a second rather than after an hour of alignment.
3. Document the parameters in `configs/template.yaml`. Every field gets a comment saying
   what it does and, where there is one, why the default is the default.
4. Write the stage contract into the subpackage's `__init__.py` docstring: inputs, outputs,
   config sections. Those docstrings are the design documentation.
5. Emit a `manifest.json` with the config hash and the versions of every external tool and
   model used. A result without provenance cannot be reproduced and should not be published.

## Conventions

- **Coordinates are 1-based and inclusive**, everywhere, matching how biologists and PDB
  files number residues. Convert at the boundary of any library that disagrees, and convert
  once.
- **Receptor numbering is canonical.** Donor coordinates are always derived through an
  explicit structural mapping, never assumed to correspond.
- **Fail loudly on partial results.** A stage that writes an incomplete output file and
  exits zero will silently poison everything downstream.
- **Record rejections, don't drop them.** Filters write their decisions to the output table
  with a reason. A filter that quietly eats most of the candidate set is a bug you can only
  see if the rejections survive.
- Structural correspondence comes from superposition, not sequence alignment. At the
  identities typical of a thermophilic donor, sequence alignment is not reliable enough to
  define a graft boundary.

## Scientific honesty

- If a third-party method's parameters are proprietary and your implementation approximates
  them, say so in the module docstring and in the output. An approximation presented as the
  published method is a misrepresentation, not a shortcut.
- Language-model scores are not stability measurements. Keep the docs clear about what the
  pipeline predicts (plausibility, direction) versus what it does not (ΔΔG, Tm).
- Defaults that encode a scientific judgement — OGT regime boundaries, score weights,
  junction rules — get a comment explaining the judgement, so a user can disagree with it
  deliberately.

## Tests

- Every config validation rule gets a test. They are cheap and they are the guard rail.
- Language-model scoring is tested against a stub implementing the `SequenceScorer`
  protocol. Do not download weights in CI.
- Pin output schemas (column names, manifest keys) in tests. Downstream steps depend on
  them, and they are the easiest thing to break during a port.
