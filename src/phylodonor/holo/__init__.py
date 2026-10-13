"""Cofactor transfer, protonation and restrained minimisation.

Stage contract
--------------
Inputs:  apo model; donor crystal PDB IDs per cofactor
Outputs: <id>_holo.pdb, <id>_qc.tsv, NN_<step>_manifest.json
Config:  holo.*
"""

from __future__ import annotations
