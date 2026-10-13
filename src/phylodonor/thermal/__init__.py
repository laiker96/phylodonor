"""Growth-temperature annotation and regime calling over clusters.

This stage decides which homologs count as thermophilic, which is what makes
one a donor candidate at all.

Why regimes are called over clusters rather than per sequence
-------------------------------------------------------------
Growth-temperature datasets join by taxon id, and that join is sometimes
wrong. Trusting each sequence's own annotation therefore puts mislabelled
mesophiles into the donor pool.

So a regime is called over a *cluster* of related sequences: a cluster needs
``min_labelled`` annotated members and ``min_purity`` agreement among them,
and one that cannot meet both is left uncalled rather than guessed. A single
bad annotation cannot carry a cluster.

This is the cheap replacement for phylogenetic clade calling. The original
pipeline got the same robustness from tree clades, which meant an
alignment/trim/tree chain -- for a benefit that never reached donor ranking,
since the donor score has no phylogenetic term. Cluster purity uses the same
mechanism, agreement among related sequences, over clusters that
``homologs`` has already computed.

What this gives up, deliberately: imputing a regime across deeper
phylogenetic distance. ``impute_unlabelled`` is false by default, because
calling a sequence thermophilic on its neighbours' annotations is a real
claim, and inheriting it silently inflates the pool with sequences nothing
is known about.

Stage contract
--------------
Inputs:  clusters.tsv and representative metadata; public OGT datasets
Outputs: sequence_ogt.tsv, cluster_regimes.tsv
Config:  thermal.*
"""

from __future__ import annotations
