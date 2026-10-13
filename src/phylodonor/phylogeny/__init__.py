"""Alignment trimming and tree inference. Opt-in; nothing depends on it.

Kept because a tree is worth having for a figure or a paper, and not because
the pipeline needs one.

Why it is not required
----------------------
Two things were checked against the original implementation. First, the
donor score is a weighted sum of per-domain structural similarity, coverage
and TM-score, with **no phylogenetic term** -- tree distances were computed
and then used only as a heatmap row-sort key, in a module nothing imported.
Second, the tree's genuine contribution was discarding noisy
growth-temperature labels by requiring agreement within a supported clade;
that is now done by cluster purity in :mod:`phylodonor.thermal`, which uses
the same mechanism far more cheaply.

Running this stage adds a trimmer and a tree inference program as
dependencies and produces nothing any other stage reads. Run it when you
want the tree itself.

Stage contract
--------------
Inputs:  alignment.fa from `homologs`
Outputs: trimmed.fa, tree.treefile, tree.contree
Config:  phylogeny.*
"""

from __future__ import annotations
