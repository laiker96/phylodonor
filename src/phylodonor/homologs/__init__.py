"""Homolog search and family classification.

Stage contract
--------------
Inputs:  target FASTA; a local sequence database
Outputs: hits.fa, hits_metadata.tsv, classification_report.tsv
Config:  search.*, classification.*
"""

from __future__ import annotations
