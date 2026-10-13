# Environment specifications

Conda/mamba specs for the optional, heavy parts of the pipeline. The core package installs
with plain `pip install -e .`; these are only needed for the stages that call external
binaries.

Ported from the source project's `config/` directory as each stage lands. Expected:

| File | For |
|---|---|
| `core.yml` | MMseqs2, HMMER, trimAl, IQ-TREE |
| `lm.yml` | PyTorch + ESM (GPU) |
| `dynamics_openmm.yml` | OpenMM + CUDA |
| `structure.yml` | gemmi, mdtraj, DSSP, US-align, Foldseek |
| `contacts.yml` | Arpeggio, GetContacts |

Licence-encumbered engines (Amber/pmemd, PyRosetta) are **not** specified here. See
[../docs/external_tools.md](../docs/external_tools.md).
