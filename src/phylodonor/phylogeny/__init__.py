"""Alignment, trimming and tree inference.

Stage contract
--------------
Inputs:  classified homolog FASTA; optional HMM profiles
Outputs: alignment.fa, column_map.tsv, tree.treefile, tree.contree
Config:  phylogeny.*, domains.hmm_profiles
"""

from __future__ import annotations
