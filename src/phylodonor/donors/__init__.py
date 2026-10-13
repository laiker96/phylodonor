"""Ranking thermophilic homologs as donors of stabilising structure.

Stage contract
--------------
Inputs:  regime-annotated tips; predicted/experimental structures; receptor holo model
Outputs: donor_scores.tsv, donor_ranking.tsv, per-domain superpositions
Config:  donors.*, domains.definitions
"""

from __future__ import annotations
