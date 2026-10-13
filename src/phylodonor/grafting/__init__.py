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

A protein language model gives a cheap, sequence-only signal for triage. The
default ranking term is:

``embedding_distance``
    How far the design's mean embedding moves toward a thermophilic reference
    set, relative to the receptor. This encodes the actual objective --
    likelihood has no idea you want thermostability -- at one forward pass per
    design. It is the default because it is the combination that has been run
    end to end here, with the surviving designs checked by structure
    prediction.

Two likelihood terms are implemented and carry no weight by default:

``delta_pll``
    Pseudo-log-likelihood of the chimera minus that of the receptor. Negative
    means the model finds the design less natural than wild type.

``junction_pll``
    The same quantity restricted to a window around each seam. Expected to be
    the most discriminating term available, since a graft can leave global
    likelihood nearly untouched while producing two implausible junctions, and
    those junctions are where designs fail experimentally.

These are better motivated than the embedding term -- masked-token
probability is what the model was trained to estimate -- but they have not
been validated in this pipeline, so the defaults do not rely on them. See
``docs/grafting.md``.

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
