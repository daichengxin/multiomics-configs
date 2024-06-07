# Cell Line Metadata Database

This repository contains scripts for creating and managing cell line metadata for the multiomics config datasets and [quantms.org](https://quantms.org), leveraging ontologies and natural language processing (NLP) to annotate cell lines in SDRF files. 

## Motivation

Cell lines are a fundamental part of biological research, and they are used in a wide range of experiments. However, cell line metadata can be inconsistent and difficult to manage. Here we are creating a [DB](cl-annotations-db.tsv) that can be used to annotate/validate proteomics SDRF for cell lines studies. These are the major sources of cell line metadata:

- [CelloSaurus](https://web.expasy.org/cellosaurus/): CelloSaurus is the main source [used in our database](cellosaurus.txt.gz). The source of the metadata can be downloaded from [cellosaurus.txt](https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt). We converted the file to a shorter version with only the fields that we are interested and the taxonomy. We use the script `pycls cellosaurus-database` to create the database.
- [Cell model passports](https://cog.sanger.ac.uk/cmp/download/model_list_20240110.csv): The cell model passports are a collection of cell lines from multiple sources. We use the file [model_list_20240110.csv](model_list_20240110.csv) to create a database extracting only the cell lines information `pycls cell-passports-to-database`.
- [EA](https://https://www.ebi.ac.uk/gxa): Expression Atlas has been curating for more than 10 years the metadata of multiple RNA experiments. We collect multiple cell lines experiments from EA in folder [ea](ea); and try to create a catalog of cell lines metadata as an extra source.
- [MONDO](https://bioportal.bioontology.org/ontologies/MONDO): The Monarch Disease Ontology (MONDO) is used to annotate the disease of the cell line.
- [BTO](https://bioportal.bioontology.org/ontologies/BTO): The BRENDA Tissue Ontology (BTO) is used to annotate an extra reference for the cell line ID. 

> **Note**: Additionally, we use other resources such as [Coriell cell line Catalog](https://www.coriell.org/), [cell bank riken](https://cell.brc.riken.jp/en/) and [atcc](https://www.atcc.org/) for manual annotation of cell lines in the database. 

## Features for every cell line

The database is created in the following path [cl-annotations-db.tsv](cl-annotations-db.tsv) and contains the following fields:

- **name**: The cell line name as defined by the curation team (ai or manual).
- **cellosaurus name**: The cell line name as annotated in Cellsaurus `ID` 
- **cellosaurus accession**: The cell line accession as annotated in Cellsaurus `AC`
- **bto cell line**: The cell line name as annotated in BTO
- **organism**: The organism of the cell line as annotated in Cellsaurus
- **organism part**: This information is not available in Cellosaurus, we use other sources to _annotate_ this field.
- **age**: The age of the cell line as annotated in Cellosaurus. If the age is not available (empty), we annotated the age from other sources such as [atcc](https://www.atcc.org/) or [Coriell cell line Catalog](https://www.coriell.org/)
- **developmental stage**: The developmental stage of the cell line as annotated in Cellosaurus; if the information is not available is inferred from the age of the cell line. 
- **sex**: Sex as provided by Cellosaurus
- **ancestry category**: The ancestry category of the cell line as annotated in Cellosaurus. If not available we use other sources. 
- **disease**: The disease is _"agreed"_ between sources CelloSaurus, EA, and ATCC.  
- **cell type**: The cell type is _"agreed"_ between sources CelloSaurus, EA, and ATCC.
- **Material**: The material is _"agreed"_ between sources CelloSaurus, EA, and ATCC.
- **synonyms**: This field is constructucted using the CelloSaurus synonyms and annotation pipelines in the library. if we find that a new synonym is available, we add it to the database.
- **curated**: This field is used to annotate if the cell line has been curated by the team, the classes are _not curated_, _ai curated_, _manual curated_.

> **Note**: The database is a tab-delimited file that can be easily read and search using pandas or GitHub table rendering. 

## Requirements

- Python 3.x
- Libraries: `pandas`, `spacy`, `click`, `owlready2`, `scikit-learn`
- Spacy model: `en_core_web_md`

## Installation

Install the required Python packages using pip:

```sh
pip install pandas spacy click owlready2 scikit-learn
python -m spacy download en_core_web_md
```

## Scripts

### pycls cell-passports-to-database

#### Summary
The `cell-passports-to-database` command reads a CSV file containing cell passport data, filters it to include only cell lines, processes and renames specific columns, and then writes the processed data to an output file in tab-separated format.

#### Example Usage

```sh 
pycls cell-passports-to-database --cell-passports path/to/cell_passports.csv --output path/to/output.tsv
```

#### Code Analysis
Inputs
- cell-passports: path to the folder containing the cell passport files
- output: path to the output file
 
#### Flow
- Read the CSV file specified by cell_passports.
- Filter the data to include only rows where model_type is "Cell Line".
- Select and rename specific columns.
- Fill missing values with "no available".
- Convert the age_at_sampling column to integers where applicable.
- Write the processed data to the specified output file in tab-separated format.
 
#### Outputs
A tab-separated file containing the processed cell passport data.