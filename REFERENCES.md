# References

The methods this pipeline builds on. Papers are cited, never redistributed — if you want
PDFs, follow the DOIs through your institution.

## Protein language models

- Rives et al. (2021) *Biological structure and function emerge from scaling unsupervised
  learning to 250 million protein sequences.* PNAS 118:e2016239118.
  doi:10.1073/pnas.2016239118 — the ESM family.
- Meier et al. (2021) *Language models enable zero-shot prediction of the effects of
  mutations on protein function.* NeurIPS 34.
  doi:10.1101/2021.07.09.450648 — the masked-marginal estimator used for `delta_pll`, and
  the ESM-1v variant-effect setup.
- Lin et al. (2023) *Evolutionary-scale prediction of atomic-level protein structure with a
  language model.* Science 379:1123. doi:10.1126/science.ade2574 — ESM-2 and ESMFold.
- Hsu et al. (2022) *Learning inverse folding from millions of predicted structures.* ICML.
  doi:10.1101/2022.04.10.487779 — ESM-IF, the optional `inverse_folding` term.
- Dauparas et al. (2022) *Robust deep learning–based protein sequence design using
  ProteinMPNN.* Science 378:49. doi:10.1126/science.add2187 — alternative inverse-folding
  scorer.

On the limits of these scores for stability specifically:

- Dieckhaus et al. (2024) *Transfer learning to leverage larger datasets for improved
  prediction of protein stability changes.* PNAS 121:e2314853121.
  doi:10.1073/pnas.2314853121 — why zero-shot likelihood is a weak ΔΔG predictor, and what
  closes the gap.

## Homology search, alignment, phylogenetics

- Steinegger & Söding (2017) *MMseqs2 enables sensitive protein sequence searching for the
  analysis of massive data sets.* Nat Biotechnol 35:1026. doi:10.1038/nbt.3988
- Eddy (2011) *Accelerated profile HMM searches.* PLoS Comput Biol 7:e1002195.
  doi:10.1371/journal.pcbi.1002195 — HMMER.
- Capella-Gutiérrez et al. (2009) *trimAl.* Bioinformatics 25:1972.
  doi:10.1093/bioinformatics/btp348
- Minh et al. (2020) *IQ-TREE 2.* Mol Biol Evol 37:1530. doi:10.1093/molbev/msaa015
- Hoang et al. (2018) *UFBoot2: improving the ultrafast bootstrap approximation.*
  Mol Biol Evol 35:518. doi:10.1093/molbev/msx281

## Growth temperature

- Sato et al. (2020) *TEMPURA: a database of growth temperatures of usable prokaryotes.*
  Microbes Environ 35:ME20074. doi:10.1264/jsme2.ME20074
- Reimer et al. (2022) *BacDive in 2022.* Nucleic Acids Res 50:D741.
  doi:10.1093/nar/gkab961

## Structure prediction and comparison

- Jumper et al. (2021) *Highly accurate protein structure prediction with AlphaFold.*
  Nature 596:583. doi:10.1038/s41586-021-03819-2
- Varadi et al. (2022) *AlphaFold Protein Structure Database.* Nucleic Acids Res 50:D439.
  doi:10.1093/nar/gkab1061
- Zhang & Skolnick (2005) *TM-align.* Nucleic Acids Res 33:2302. doi:10.1093/nar/gki524
- Zhang et al. (2022) *US-align: universal structure alignment.* Nat Methods 19:1109.
  doi:10.1038/s41592-022-01585-1
- van Kempen et al. (2024) *Fast and accurate protein structure search with Foldseek.*
  Nat Biotechnol 42:243. doi:10.1038/s41587-023-01773-0
- Kabsch & Sander (1983) *Dictionary of protein secondary structure.* Biopolymers 22:2577.
  doi:10.1002/bip.360221211 — DSSP, used for junction placement.
- Jubb et al. (2017) *Arpeggio: a web server for calculating and visualising interatomic
  interactions in protein structures.* J Mol Biol 429:365. doi:10.1016/j.jmb.2016.12.004

## Molecular dynamics

- Eastman et al. (2017) *OpenMM 7.* PLoS Comput Biol 13:e1005659.
  doi:10.1371/journal.pcbi.1005659
- Case et al. *AmberTools.* J Chem Inf Model 63:6183 (2023).
  doi:10.1021/acs.jcim.3c00644
- Maier et al. (2015) *ff14SB.* J Chem Theory Comput 11:3696. doi:10.1021/acs.jctc.5b00255

## Solubility and aggregation

- Sormanni et al. (2015) *The CamSol method of rational design of protein mutants with
  enhanced solubility.* J Mol Biol 427:478. doi:10.1016/j.jmb.2014.09.026
- Kuriata et al. (2019) *Aggrescan3D standalone package.* Bioinformatics 35:3834.
  doi:10.1093/bioinformatics/btz143

## Thermostability engineering by consensus and ancestral approaches

Context for why donor grafting works at all, and the main alternatives to it:

- Lehmann & Wyss (2001) *Engineering proteins for thermostability: the use of sequence
  alignments versus rational design and directed evolution.* Curr Opin Biotechnol 12:371.
  doi:10.1016/S0958-1669(00)00229-9
- Perl et al. (2000) *Two exposed amino acid residues confer thermostability on a cold
  shock protein.* Nat Struct Biol 7:380. doi:10.1038/75151
- Risso et al. (2013) *Hyperstability and substrate promiscuity in laboratory resurrections
  of Precambrian β-lactamases.* JACS 135:2899. doi:10.1021/ja311630a
- Wijma et al. (2013) *Computationally designed libraries for rapid enzyme
  stabilization.* Protein Eng Des Sel 26:709. doi:10.1093/protein/gzt057
