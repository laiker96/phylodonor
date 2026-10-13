"""Homolog search, subfamily filtering, clustering and alignment.

Four jobs in one stage, because they are a single pass over the same data and
the later ones are cheap:

**search** finds family members in a sequence database.

**classify** filters to the right subfamily. Not optional: a donor from a
different subfamily may share a fold name and not a mechanism, and mixing
subfamilies also corrupts the alignment.

**cluster** keeps one representative per cluster. This stops an
over-sequenced genus from dominating everything downstream, and the clusters
become the grouping unit that ``thermal`` calls growth-temperature regimes
over -- which is how label noise gets filtered without building a tree.

**align** produces a multiple sequence alignment of the representatives,
needed for the ``donor_conservation`` term in ``grafting``. No trimming and
no tree; see ``docs/roadmap.md``.

Stage contract
--------------
Inputs:  target FASTA; a local sequence database; optional HMM profiles
Outputs: hits.fa, hits_metadata.tsv, classification_report.tsv,
         clusters.tsv, representatives.fa, alignment.fa, column_map.tsv
Config:  search.*, classification.*, cluster.*, align.*
"""

from __future__ import annotations
