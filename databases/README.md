## Databases

This directory contains the databases used in the project. The databases are stored in fasta files with extensions fasta or fa.

### Contaminats databases

- [contaminants.fasta](contaminants.fasta): A database of common contaminants in proteomics experiments. This database is used to filter out common contaminants from the search results. It is the merge of crap-202105.fasta and contaminants-mq-202105.fasta from the [MaxQuant](https://www.maxquant.org/) software.

- [contaminants-202105-uniprot.fasta](contaminants-202105-uniprot.fasta): Based on the contaminants.fasta file, this database contains only the contaminants that have a UniProt accession number.

- [contaminants-202105-uniprot-description.fasta](contaminants-202105-uniprot-description.fasta): Based on the contaminants-202105-uniprot.fasta file, this database contains the contaminants with the UniProt accession number and the protein description.

