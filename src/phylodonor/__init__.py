"""phylodonor -- phylogeny-guided donor grafting for protein thermostability.

Stages, in pipeline order. Required:

    homologs   -> search, filter to subfamily, cluster, align
    thermal    -> annotate growth temperature, call regimes over clusters
    donors     -> rank thermophilic homologs as donors
    holo       -> build cofactor-bound structures, annotate contacts
    grafting   -> enumerate chimeric designs, triage them with a protein LM

Optional -- nothing above reads their output:

    scoring    -> conservation, solubility, aggregation predictors
    viz        -> logos, plots
    phylogeny  -> alignment trimming and tree inference
    dynamics   -> molecular dynamics

A tree contributes nothing to donor ranking: the donor score has no
phylogenetic term. Its one real contribution -- discarding noisy
growth-temperature labels -- is done by cluster purity in ``thermal``
instead, which reuses the clustering ``homologs`` already performs. Likewise
no stage needs molecular dynamics. See ``docs/roadmap.md``.

Each stage reads its parameters from a single YAML config; see
``phylodonor.config`` and ``configs/template.yaml``.
"""

from __future__ import annotations

__version__ = "0.1.0.dev0"

from phylodonor.config import Config, ConfigError, load  # noqa: F401

__all__ = ["Config", "ConfigError", "load", "__version__"]
