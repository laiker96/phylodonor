"""Cofactor-bound model building: transfer, minimisation and contact annotation.

Why this stage exists
---------------------
Structure predictors give you an apo model. For an enzyme whose function
depends on a cofactor that is the wrong starting point: the cofactor pocket
has nothing in it, so minimisation and dynamics collapse it, and every
contact-based analysis downstream measures a site that does not exist in the
real protein.

Experimental holo structures of close homologs do exist, so the cofactor's
coordinates can be borrowed from one. That is what this stage does, plus the
bookkeeping needed to make the result usable by a force field and by the
later stages.

Pipeline
--------
1. **transfer** -- normalise the receptor to a single protein chain; for each
   cofactor, superpose a donor crystal by sequence-guided CA alignment and
   copy the cofactor's heavy atoms into the receptor frame.
2. **minimize** -- parameterise and relax the assembled model with the
   backbone restrained, so the cofactor settles into the pocket without the
   fold drifting.
3. **contacts** -- annotate which receptor residues touch which ligand, with
   clash detection.

Annotation ligands -- nucleic acids, substrate analogues, inhibitors present
in the donor crystal -- are transferred into a *separate* model used only for
contact calling. They are evidence about where the substrate binds, not part
of the enzyme, and including them in the minimised system would be a claim
nobody made.

Stage contract
--------------
Inputs
    - ``target.structure``, or a predicted model if that is null
    - ``holo.cofactors``: resname plus a donor PDB id for each
    - ``holo.annotation_ligands``: optional, for contact annotation only
Outputs
    - ``<model_id>_holo.pdb``          protein + cofactors, for minimisation
    - ``<model_id>_annotation.pdb``    protein + annotation ligands
    - ``<model_id>_minimized.pdb``     plus standard-resname and no-H variants
    - ``<model_id>_contacts.tsv``      per-residue contact annotation
    - ``manifest.json``                donors, provenance, fit quality, versions
Config
    ``holo.*``, ``target.structure``

Downstream, ``grafting`` reads the contact table to build its no-go list, so
residue numbering here is the receptor's throughout and is never renumbered.

A note on reading the fit quality
---------------------------------
Every superposition's RMSD and pair count go into the manifest, and they are
worth looking at rather than skipping. A global fit over few pairs, or with a
high RMSD, means the donor is not as close a homolog as assumed, and the
cofactor placement inherits that error. The stage refuses outright below
``min_fit_pairs`` -- a transform from a handful of CA atoms is arbitrary, and
returning it silently would produce a model that looks fine and is wrong.

Modules
-------
``structure``
    Parsing, normalisation, residue-name canonicalisation, PDB record I/O.
``superpose``
    Sequence-guided CA superposition with outlier rejection, and the
    pocket-restricted variant.
``fetch``
    Donor crystal retrieval from the RCSB PDB, cached per run.
``transfer``
    Cofactor and annotation-ligand transfer. The scientific core.
``contacts``
    Per-residue contact annotation with clash and hydrogen-bond heuristics.
"""

from __future__ import annotations

__all__ = ["contacts", "fetch", "structure", "superpose", "transfer"]
