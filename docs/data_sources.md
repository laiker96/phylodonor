# Reference data

Nothing in this list is vendored or downloaded automatically. You obtain each dataset and
point the config at your local copy. That is deliberate: these are large, they are
third-party, and several have redistribution terms that a toolkit repository should not
quietly assume.

## Sequence databases

For `search.database`. Any MMseqs2-indexable database works.

| Database | Notes |
|---|---|
| UniRef90 / UniRef50 | The usual choice. UniRef50 is much faster and adequate for finding distant thermophilic relatives; UniRef90 gives better clade resolution. |
| UniProt (SwissProt + TrEMBL) | Use when you need full annotation on every hit. |
| MGnify / JGI metagenomes | More thermophilic diversity, but organism-level growth-temperature annotation is usually unavailable, which weakens the `thermal` stage. |

Build an index once:

```bash
mmseqs createdb uniref90.fasta /data/db/uniref90
mmseqs createindex /data/db/uniref90 /tmp
```

Then set `search.database: /data/db/uniref90`.

## Growth-temperature datasets

For `thermal.sources`. These map an organism or NCBI taxon ID to an optimal growth
temperature, which is what turns a tree into a thermal hypothesis. Use more than one —
coverage is patchy and they disagree.

| Source | Coverage | Citation required |
|---|---|---|
| **TEMPURA** | ~8.6 k prokaryotes, curated optimal growth temperatures | Sato et al. 2020 |
| **OGT predictions from genomes** | wider, lower confidence; predicted rather than measured | cite the specific dataset |
| **BacDive** | large curated phenotype database; API access, registration required | cite BacDive and the retrieval date |

Two cautions that matter more than they sound:

- **Measured and predicted OGT must not be mixed silently.** A clade called thermophilic on
  predicted values is a weaker claim than one called on measured values, and the whole donor
  ranking rests on that call. Keep the provenance column.
- **Growth temperature is not protein melting temperature.** A thermophile's proteins are
  usually more stable, but the relationship is loose and per-protein. OGT is a prior for
  *where to look*, not a measurement of the donor.

## Structures

For the receptor and for donors. Predicted structures are fine and usually necessary —
most homologs have no experimental structure.

| Source | Use |
|---|---|
| **AlphaFold DB** | bulk donor models, by UniProt accession. Check pLDDT per region before trusting a junction placement. |
| **RCSB PDB** | experimental structures; required as cofactor donors in the `holo` stage. |
| **Local prediction** | AlphaFold, Boltz, ESMFold for sequences absent from AFDB, and for designs. |

A low-pLDDT loop in a donor model is not a reason to reject the donor, but it *is* a reason
to reject junctions placed inside that loop. The `holo` stage emits per-residue confidence
so `grafting` can use it.

## Cofactor donor crystals

For `holo.cofactors`. Each cofactor needs an experimental structure to transfer coordinates
from — a high-resolution holo structure of a close homolog, chosen for resolution and for
having the cofactor fully ordered. Predicted structures generally place cofactors poorly or
not at all, which is why this stage exists.

## Citation

If you publish, cite every dataset you used, with the version or access date. Most of these
are maintained by small groups for whom citations are the entire return on the work.
