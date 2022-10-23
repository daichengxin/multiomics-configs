# Detection of non-canonical peptides using custom databases

We have implemented the pypgatk package and the pgdb workflow to create proteogenomics databases based on ENSEMBL resources. The tools allow the generation of protein sequences from novel protein-coding transcripts by performing a three-frame translation of pseudogenes, lncRNAs, and other non-canonical transcripts, such as those produced by alternative splicing events. It also includes exonic out-of-frame translation from otherwise canonical protein-coding mRNAs. Moreover, the tool enables the generation of variant protein sequences from multiple sources of genomic variants including COSMIC, cBioportal, gnomAD, and mutations detected from sequencing of patient samples.

The preprint of the original manuscript can be found here [Generation of ENSEMBL-based proteogenomics databases boosts the identification of non-canonical peptides](https://www.biorxiv.org/content/10.1101/2021.06.08.447496v1).

List of Datasets analyzed:

- [PXD005946](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD005946) [61 Samples] - Global Proteome Analysis of the NCI-60 Cell Line Panel, part 3
- [PXD019263](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD005946) [2 Samples] - HepG2 whole cell lysate LC-MS and 2D-LC-MS
- [PXD004452](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD004452) [4 Samples]- HeLa proteome of 12,250 protein-coding genes
- [PXD014145](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD014145) [1 Samples]- Depletion of histone methyltransferase KMT9 inhibits lung cancer cell proliferation by inducing non-apoptotic cell death.
- [PXD020467](To be annotated)
- [PXD020454](To be annotated)
- 