"""phylodonor -- phylogeny-guided donor grafting for protein thermostability.

Stages, in pipeline order:

    homologs   -> find and classify family members
    phylogeny  -> align, trim, build a tree
    thermal    -> annotate growth temperature, call thermal-regime clades
    donors     -> rank thermophilic homologs as donors
    holo       -> build cofactor-bound, minimised structures
    grafting   -> enumerate chimeric designs, rank them with a protein LM
    dynamics   -> simulate survivors
    scoring    -> sequence-level predictors (conservation, solubility, aggregation)
    viz        -> trees, logos, plots

Each stage reads its parameters from a single YAML config; see
``phylodonor.config`` and ``configs/template.yaml``.
"""

from __future__ import annotations

__version__ = "0.1.0.dev0"

from phylodonor.config import Config, ConfigError, load  # noqa: F401

__all__ = ["Config", "ConfigError", "load", "__version__"]
