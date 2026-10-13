"""Sequence-level predictors applied to receptor and designs.

Stage contract
--------------
Inputs:  sequences; alignment (for conservation)
Outputs: conservation.tsv, solubility.tsv, aggregation.tsv
Config:  scoring.*
"""

from __future__ import annotations
