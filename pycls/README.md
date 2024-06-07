# Cell Line Metadata Database

This repository contains scripts for creating and managing cell line metadata for the multiomics config datasets and [quantms.org](https://quantms.org), leveraging ontologies and natural language processing (NLP) to annotate cell lines in SDRF files. 

## Motivation

Cell lines are a fundamental part of biological research, and they are used in a wide range of experiments. However, cell line metadata can be inconsistent and difficult to manage. Here we are creating a [DB](cl-annotations-db.tsv) that can be used to annotate/validate proteomics SDRF for cell lines studies. These are the major sources of cell line metadata:

- [CelloSaurus](https://web.expasy.org/cellosaurus/): CelloSaurus is the main source used in our database. The source of the metadata is the following file [cellosaurus.txt](https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt)
- [Cell model passports](https://cog.sanger.ac.uk/cmp/download/model_list_20240110.csv): The cell model passports are a collection of cell lines from multiple sources. We use the cell model passports to annotate cell lines that are not available in CelloSaurus.
- [BTO](https://bioportal.bioontology.org/ontologies/BTO): The BRENDA Tissue Ontology (BTO) is used to annotate an extra reference for the cell line ID. 
- [EA](https://https://www.ebi.ac.uk/gxa): Expression Atlas has been curating for more than 10 years the metadata of multiple RNA experiments. We collect multiple cell lines experiments from EA in folder [ea](ea); and try to create a catalog of cell lines metadata as an extra source.
- [MONDO](https://bioportal.bioontology.org/ontologies/MONDO): The Monarch Disease Ontology (MONDO) is used to annotate the disease of the cell line.

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

### 1. Create Cell Line Metadata Database

This script creates a cell line metadata database using the CelloSaurus database, BTO and MONDO ontologies, and SDRF files.

#### Usage

```sh
python script_name.py cl-database --cellosaurus <cellosaurus_file> --bto <bto_file> --mondo <mondo_file> --sdrf_path <sdrf_folder> --database <database_file> --unknown <unknown_file>
```

#### Options

- `--cellosaurus`: Path to the gzipped CelloSaurus database file.
- `--bto`: Path to the BTO ontology file.
- `--mondo`: Path to the MONDO ontology file.
- `--sdrf_path`: Folder containing SDRF files.
- `--database`: Path to the existing database file with cell line metadata.
- `--unknown`: Output file for unknown cell lines.

### 2. NLP-Based Cell Line Recommendation

This script uses NLP to find cell lines in the CelloSaurus database when the names do not match exactly.

#### Usage

```sh
python script_name.py nlp-recommendation --cellosaurus <cellosaurus_file> --unknown <unknown_file> --output <output_file>
```

#### Options

- `--cellosaurus`: Path to the gzipped CelloSaurus database file.
- `--unknown`: File with unknown cell lines.
- `--output`: Output file with the recommendations.

### 3. Create Expression Atlas Database

This script creates a database from large Expression Atlas files.

#### Usage

```sh
python script_name.py ea-database --ea-folder <ea_folder> --output <output_file>
```

#### Options

- `--ea-folder`: Folder containing Expression Atlas files.
- `--output`: Output file with the database.

## Functions

### Parsing Functions

- **parse_cellosaurus_file(file_path)**: Parses the CelloSaurus file and returns a list of dictionaries with the parsed content.
- **read_obo_file(file_path)**: Reads an OBO file and returns a list of dictionaries with the parsed content.
- **read_cell_line_database(database)**: Reads the cell line database and returns the data as a list of dictionaries.
- **write_database(current_cl_database, comments, database)**: Writes the current cell line database list to the database file.
- **get_cell_line_code(sdrf_file)**: Extracts cell line codes from SDRF files.
- **parse_cellosaurus_taxonomy(organism_text)**: Parses the organism text from the CelloSaurus database.
- **is_age_in_text(age_text)**: Checks if the age field contains a number indicating an age rather than a developmental stage.

### Context Functions

- **modo_dict_to_context(obo_list)**: Converts OBO dictionaries to context strings for NLP processing.
- **cellosaurus_dict_to_context(cellosaurus_list)**: Converts CelloSaurus dictionaries to context strings for NLP processing.

### Similarity Functions

- **calculate_similarity(cellosaurus_text, synonyms)**: Calculates the cosine similarity between the cellosaurus text and synonyms.
- **map_celllines(cellosaurus_text, context)**: Maps cell lines using the context and NLP model.

### Entry Creation Functions

- **create_new_entry(old_cl, bto, cellosaurus_list, modo_context)**: Creates a new entry for a cell line not found in the database.
- **find_cell_line(old_cl, current_cl_database)**: Finds a given cell line in the standardized cell line database.
- **get_cell_line_bto(bto_code, bto_list)**: Retrieves the BTO cell line name from the BTO code.
- **is_in_synonyms(old_cl, cellosaurus)**: Checks if a cell line is in the synonyms list.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

Feel free to replace `script_name.py` with the actual name of your script file. Adjust the sections as necessary to fit your specific use case and additional details.