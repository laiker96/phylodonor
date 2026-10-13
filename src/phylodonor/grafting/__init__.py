"""Donor-to-receptor grafting, ranked by a protein language model.

This is the design stage. Everything upstream of it is reconnaissance: find
homologs, place them on a tree, annotate growth temperature, pick donors that
live hotter than the receptor. This stage decides *what to actually build*.

The receptor is the target protein. A donor is a thermophilic homolog. A
**graft** replaces one or more structurally equivalent segments of the
receptor with the donor's version of them.

Why a language model
--------------------
The enumeration step is generous on purpose -- a few structurally plausible
regions across a handful of donors produces hundreds to low thousands of
candidate chimeras. Folding and simulating all of them is not affordable, and
structural filters alone do not distinguish a graft that reads like a real
protein from one that merely fits geometrically.

A protein language model gives a cheap, sequence-only prior for "does this
look like something evolution would produce", which is a useful proxy for
foldability and expressibility. Two terms do most of the work:

``delta_pll``
    Pseudo-log-likelihood of the chimera minus that of the unmodified
    receptor. Negative means the model finds the design less natural than
    wild type. Small negative values are normal and acceptable; large drops
    almost always indicate a broken graft.

``junction_pll``
    The same quantity restricted to a window around each seam. This is the
    discriminating term in practice. A graft can leave global likelihood
    nearly untouched while producing two implausible junctions, and those
    junctions are exactly where designs fail experimentally.

Optionally, an inverse-folding model (ESM-IF, ProteinMPNN) scores the chimeric
sequence *conditioned on the receptor backbone*, which asks a sharper
question: not "is this a plausible protein" but "is this a plausible protein
for this fold". It costs a structure prediction per design, so it is off by
default and worth enabling once the candidate set is down to tens.

What the language model does not tell you
-----------------------------------------
Likelihood is not stability. The model has no notion of your assay, your
expression host, or your cofactor. It reliably rejects nonsense and it ranks
survivors sensibly, but the ordering within the top handful is not
trustworthy enough to skip simulation or experiment. Treat this stage as a
filter that takes thousands of designs to tens, not as an oracle that picks
the winner.

Stage contract
--------------
Inputs
    - receptor sequence and structure (``target`` in the config)
    - ranked donors with per-domain structural alignments (``donors`` stage)
    - cofactor contacts and secondary structure for the receptor (``holo``)
Outputs
    - ``designs.tsv``   one row per candidate, all score terms, ranked
    - ``designs.fa``    sequences of the surviving designs
    - ``manifest.json`` config hash, model versions, tool versions

Modules
-------
``models``
    Dataclasses shared across the stage.
``candidates``
    Candidate generation and the structural/positional constraints.
``score``
    Language-model scoring of candidates.
``rank``
    Composite scoring, filtering, and output writing.
"""

from __future__ import annotations

__all__ = ["models", "candidates", "score", "rank"]
