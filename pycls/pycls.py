import gzip
from typing import Union

import click
import pandas as pd
import glob

import spacy
import numpy as np
from owlready2 import *
from sklearn.metrics.pairwise import cosine_similarity

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

nlp = spacy.load("en_core_web_md")  # Load the spacy model


def get_cell_line_code(sdrf_file):
    sdrf = pd.read_csv(sdrf_file, sep="\t")
    cl_list = sdrf["characteristics[cell line]"].unique().tolist()
    return cl_list


def parse_cellosaurus_file(file_path):
    """
    Parse the CelloSaurus file and return a list of dictionaries with the parsed content
    :param file_path: CelloSaurus file path
    :return: List of CelloSaurus dictionaries
    """

    def parse_entry(entry):
        data = {}
        if "ID   FU-OV-1" in entry:
            print(entry)
        lines = entry.strip().split("\n")
        for line in lines:
            if line.startswith("ID"):
                data["name"] = line.split("ID ")[1].strip()
            elif line.startswith("AC"):
                data["accession"] = line.split("AC ")[1].strip()
            elif line.startswith("SY"):
                data["synonyms"] = line.split("SY ")[1].strip().split("; ")
            elif line.startswith("DR   BTO"):
                data["bto"] = line.split("; ")[1]
            elif line.startswith("DR   EFO"):
                data["efo"] = line.split("; ")[1]
            elif line.startswith("OX"):
                data["organism"] = line.split("OX ")[1].strip()
            elif line.startswith("SX"):
                data["sex"] = line.split()[1]
            elif line.startswith("AG"):
                data["age"] = line.split("AG ")[1].strip()
            elif line.startswith("CA"):
                data["cell_type"] = line.split()[1]
            elif line.startswith("CC") and "Population" in line:
                data["ancestry category"] = line.split(": ")[1].strip()
            elif line.startswith("DI"):
                if "NCIt" in line:
                    # Regular expression to match all disease annotations and capture the disease name
                    pattern = r"NCIt;\s*C\d+;\s*([^;]+)"
                    match = re.search(pattern, line)
                    if match:
                        data["disease"] = match.group(1)

        return data

    # Read the file and split into entries, the file is gzipped
    with gzip.open(file_path, "r") as file:
        content = file.read()

    # Split the content by entries
    entries = content.split("//\n")

    # Parse each entry
    parsed_data = [parse_entry(entry) for entry in entries if entry.strip()]
    # remove empty entries
    parsed_data = [entry for entry in parsed_data if entry]
    return parsed_data


def read_obo_file(file_path):
    """
    Read an obo file and return a list of dictionaries with the parsed content
    :param file_path: OBO file path
    :return: List of OBO dictionaries
    """
    with open(file_path, "r") as file:
        content = file.read()

    # Split the content by entries
    entries = content.split("\n\n")

    def parse_obo_term(entry: str) -> dict:
        obo_dict = {}
        lines = entry.strip().split("\n")
        for line in lines:
            if line.startswith("id:"):
                # extract the id from pattern id: BTO:0000001 using regex
                obo_dict["id"] = line.split("id: ")[1].strip()
            elif line.startswith("name:"):
                obo_dict["name"] = line.split("name: ")[1].strip()
            elif line.startswith("def:"):
                obo_dict["def"] = line.split("def: ")[1].strip()
            elif line.startswith("is_a:"):
                if "is_a" not in obo_dict:
                    obo_dict["is_a"] = []
                obo_dict["is_a"].append(line.split("is_a: ")[1].strip())
            elif line.startswith("synonym:"):
                if "synonyms" not in obo_dict:
                    obo_dict["synonyms"] = []
                obo_dict["synonyms"].append(line.split("synonym: ")[1].strip())
        return obo_dict

    obo_list = [parse_obo_term(entry) for entry in entries if entry.strip()]

    # remove empty entries
    obo_list = [entry for entry in obo_list if entry]

    return obo_list


def modo_dict_to_context(obo_list: list) -> str:
    context = ""
    for entry in obo_list:
        if "obsolete" not in entry["name"]:
            context += f"{entry['id']}: {entry['name']}\n"
    return context


def calculate_similarity(cellosaurus_text, synonyms):
    query_vec = nlp(" ".join(cellosaurus_text)).vector
    synonyms_vec = [nlp(" ".join(synonym)).vector for synonym in synonyms]
    similarities = cosine_similarity([query_vec], synonyms_vec)
    return max(similarities[0])


def map_celllines(cellosaurus_text: str, context: list):
    # Use the LLM to find the correct MONDO term
    max_similarity = 0
    closest_match = None
    for entry in context:
        synonyms = entry.split(";")
        similarity = calculate_similarity(cellosaurus_text, synonyms)
        if similarity > max_similarity:
            closest_match = entry
            max_similarity = similarity
    return closest_match


def read_cell_line_database(database):
    """
    The database is a dictionary with the following structure. First comments are ignored.
    Each block of information is devided by //. The information for each cell lines is:

    cellosaurus name
    bto cell line
    organism
    organism part
    age
    developmental stage
    sex
    ancestry category
    disease
    cell type
    Material
    synonyms

    The cell line is a BTO name. and the synonyms are separated by ; and the disease is a MONDO disease term.
    :param database:
    :return:
    """

    with open(database, "r") as file:
        content = file.read()

    # get comments and ignore them
    comments = [line for line in content.split("\n") if line.startswith("#")]

    # Split the content by entries
    entries = content.split("//")

    def parse_database_entry(entry: str) -> dict:
        data = {}
        lines = entry.strip().split("\n")
        for line in lines:
            if line.startswith("cellosaurus name"):
                data["cellosaurus name"] = line.split(": ")[1].strip()
            elif line.startswith("bto cell line"):
                data["cell line"] = line.split(": ")[1].strip()
            elif line.startswith("organism"):
                data["organism"] = line.split(": ")[1].strip()
            elif line.startswith("organism part"):
                data["organism part"] = line.split(": ")[1].strip()
            elif line.startswith("age"):
                data["age"] = line.split(": ")[1].strip()
            elif line.startswith("developmental stage"):
                data["developmental stage"] = line.split(": ")[1].strip()
            elif line.startswith("sex"):
                data["sex"] = line.split(": ")[1].strip()
            elif line.startswith("ancestry category"):
                data["ancestry category"] = line.split(": ")[1].strip()
            elif line.startswith("disease"):
                data["disease"] = line.split(": ")[1].strip()
            elif line.startswith("cell type"):
                data["cell type"] = line.split(": ")[1].strip()
            elif line.startswith("Material"):
                data["Material"] = line.split(": ")[1].strip()
            elif line.startswith("synonyms"):
                data["synonyms"] = line.split("synonyms: ")[1].strip().split("; ")
        return data

    database_list = [parse_database_entry(entry) for entry in entries if entry.strip()]
    # remove empty entries
    database_list = [entry for entry in database_list if entry]
    return database_list, comments


def find_cell_line(old_cl: str, current_cl_database: list) -> Union[dict, None]:
    """
    Find a given cell line annotated in an SDRF in an standarized cell line database
    :param old_cl: Code (e.g. HELA, A549, etc.)
    :param current_cl_database: Database of all cell lines for the multiomics configuration
    :return:
    """

    # Normalize the cell line name to lower case and remove spaces
    old_cl = old_cl.lower().strip()

    for entry in current_cl_database:
        if "cell line" in entry:
            if entry["cell line"].lower().strip() == old_cl:
                return entry
        if "cellosaurus name" in entry:
            if entry["cellosaurus name"].lower().strip() == old_cl:
                return entry
        for synonym in entry["synonyms"]:
            if synonym.lower().strip() == old_cl:
                return entry
    return None


def is_in_synonyms(old_cl, cellosaurus) -> bool:
    """
    Check if a cell line is in the synonyms list
    :param old_cl: Old cell line code
    :param sysnonyms: Synonyms list
    :return:
    """
    if "synonyms" not in cellosaurus:
        return False
    for synonym in cellosaurus["synonyms"]:
        if synonym.lower().strip() == old_cl.lower().strip():
            return True
    return False


def get_cell_line_bto(bto_code: str, bto_list: list) -> Union[None, str]:
    for entry in bto_list:
        if entry["id"] == bto_code:
            return entry["name"]
    return None


def parse_cellosaurus_taxonomy(organism_text: str):
    """
    Parse the organism text from the cellosaurus database
    :param organism_text: organism text from the cellosaurus database
    :return: return the species name and the taxonomy id
    """
    pattern = r"NCBI_TaxID=(\d+); ! ([\w\s]+) \(([\w\s]+)\)"
    match = re.search(pattern, organism_text)

    if match:
        ncbi_tax_id = match.group(1)
        species_name = match.group(2)
        return species_name, ncbi_tax_id
    return None, None


def is_age_in_text(age_text: str) -> bool:
    """
    Check if the age field contains a number is an Age if not is a developmental stage.
    :param age_text: String text
    :return: True if the age is a number, False otherwise
    """

    return any(char.isdigit() for char in age_text)


def create_new_entry(old_cl, bto, cellosaurus_list, modo_context):
    """
    Create a new entry for a cell line not found in the database
    :param old_cl:
    :param bto:
    :param cellosaurus_list:
    :param modo_context:
    :return:
    """
    old_cl = old_cl.lower().strip()
    for cellosaurus in cellosaurus_list:
        if "name" in cellosaurus and (
            cellosaurus["name"].lower().strip() == old_cl
            or is_in_synonyms(old_cl, cellosaurus)
        ):
            entry = {}
            entry["cellosaurus name"] = cellosaurus["name"]
            if "bto" in cellosaurus:
                bto_name = get_cell_line_bto(cellosaurus["bto"], bto)
                if bto_name is not None:
                    entry["bto cell line"] = bto_name
                else:
                    raise ValueError(
                        f"Cell line {old_cl} not found in the BTO ontology"
                    )
            if "organism" in cellosaurus:
                scientific_name, tax = parse_cellosaurus_taxonomy(
                    cellosaurus["organism"]
                )
                if scientific_name is not None:
                    entry["organism"] = scientific_name
                    entry["taxonomy"] = tax
            entry["organism part"] = cellosaurus["cell_type"]
            if "age" in cellosaurus:
                if is_age_in_text(cellosaurus["age"]):
                    entry["age"] = cellosaurus["age"]
                    entry["developmental stage"] = None
                else:
                    entry["developmental stage"] = cellosaurus["age"]
                    entry["age"] = None
            if "sex" in cellosaurus:
                entry["sex"] = cellosaurus["sex"]
            else:
                entry["sex"] = None
            if "ancestry category" in cellosaurus:
                entry["ancestry category"] = cellosaurus["ancestry category"]
            else:
                entry["ancestry category"] = None
            if "disease" in cellosaurus:
                entry["disease"] = cellosaurus["disease"]
            else:
                entry["disease"] = None
            # entry["bto"] = bto
            # entry["modo"] = modo_context
            if "synonyms" in cellosaurus:
                entry["synonyms"] = cellosaurus["synonyms"]
            else:
                entry["synonyms"] = []
            return entry
    return None


def write_database(current_cl_database: list, comments: list, database: str):
    """
    Write the database objects to the database file
    :param current_cl_database: current cell line database list
    :param comments: comments in the database file
    :param database: database file path
    :return:
    """

    with open(database, "w") as file:
        for comment in comments:
            file.write(comment + "\n")
        for entry in current_cl_database:
            file.write("//\n")
            for key, value in entry.items():
                if key == "synonyms":
                    if value:
                        value = "; ".join(value)
                    else:
                        value = None
                file.write(f"{key}: {value}\n")


def cellosaurus_dict_to_context(cellosaurus_list: list) -> list:
    context = []
    for entry in cellosaurus_list:
        if "synonyms" in entry:
            synonyms = "; ".join(entry["synonyms"])
        if "name" in entry:
            context.append(f"{synonyms}; {entry['name']}")
    context = [preprocess_text(entry) for entry in context]
    print("Cellosaurus context created -- ", len(context))
    return context


def preprocess_text(text):
    # Tokenize and preprocess text
    return text.lower().strip()


@click.command(
    "cl-database",
    short_help="Create a cell lines metadata database for annotating cell lines sdrfs",
)
@click.option(
    "--cellosaurus",
    help="CelloSaurus database file, the file is gzipped",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "--bto", help="BTO ontology file", required=True, type=click.Path(exists=True)
)
@click.option(
    "--mondo", help="Mondo ontology file", required=True, type=click.Path(exists=True)
)
@click.option(
    "--sdrf_path", help="SDRF folder with all existing SDRF files", required=True
)
@click.option(
    "--database", help="Existing database file with cell lines metadata", required=True
)
@click.option(
    "--unknown", help="Output for unknown cell lines in cellosaurus", required=True
)
def cl_database(
    cellosaurus: str, bto: str, mondo: str, sdrf_path: str, database: str, unknown: str
) -> None:
    """
    The following function creates a vector database using LLMs using the CelloSaurus database, BTO and EFO ontologies
    :param cellosaurus: CelloSaurus database file
    :param bto: BTO ontology file
    :param efo: EFO ontology file
    :param output: Output file with a vector database constructed using LLMs
    :return:
    """

    # Read the current cell line database
    current_cl_database, comments = read_cell_line_database(database)

    # Parse the CelloSaurus file
    cellosaurus_list = parse_cellosaurus_file(cellosaurus)

    # Read the BTO and EFO from owl files
    bto = read_obo_file(bto)

    # Read the EFO from owl files
    modo = read_obo_file(mondo)
    modo_context = modo_dict_to_context(modo)

    sdrf_files = glob.glob(sdrf_path + "/**/*.tsv", recursive=True)

    cl_list = []
    for sdrf_file in sdrf_files:
        cl_list += get_cell_line_code(sdrf_file)

    cl_list = list(set(cl_list))
    print("Number of cell lines in the SDRF files: ", len(cl_list))

    # Check the database if the cell line is not already, create a new entry
    # If the cell is ready to get the information from the database.

    non_found_cl = []
    for old_cl in cl_list:
        cl_db = find_cell_line(old_cl, current_cl_database)
        if not cl_db:
            # print(f"Cell line {old_cl} not found in the database - attend to create one programmatically")
            cl_new_entry = create_new_entry(
                old_cl, bto, cellosaurus_list, modo_context
            )  # Create a new entry
            if cl_new_entry is not None:
                if current_cl_database is None:
                    current_cl_database = []
                current_cl_database.append(cl_new_entry)
            else:
                print(f"{old_cl}")
                non_found_cl.append(old_cl)
        else:
            print("Cell line found in the database: ", cl_db["cellosaurus name"])

    write_database(current_cl_database, comments, database)

    # Write the unknown cell lines to a file
    with open(unknown, "w") as file:
        for cl in non_found_cl:
            file.write(cl + "\n")


@click.command(
    "nlp-recommendation",
    short_help="Looks for cell lines in cellosaurus using NLP because the name do not match the cellosaurus name",
)
@click.option(
    "--cellosaurus",
    help="CelloSaurus database file",
    required=True,
    type=click.Path(exists=True),
)
@click.option("--unknown", help="unknown cell lines", required=True)
@click.option("--output", help="File with the recomendations", required=True)
def nlp_recommendation(cellosaurus: str, unknown: str, output: str) -> None:
    """
    The following function creates a vector database using LLMs using the CelloSaurus database, BTO and EFO ontologies
    :param cellosaurus: CelloSaurus database file
    :param bto: BTO ontology file
    :param efo: EFO ontology file
    :param output: Output file with a vector database constructed using LLMs
    :return:
    """

    # Parse the CelloSaurus file
    cellosaurus_list = parse_cellosaurus_file(cellosaurus)

    cellosaurus_context = cellosaurus_dict_to_context(cellosaurus_list)

    with open(unknown, "r") as file:
        content = file.read()

    unknown_cl = content.split("\n")

    with open(output, "w") as file:
        for cl in unknown_cl:
            search_term = preprocess_text(cl)
            llm_term = map_celllines(search_term, cellosaurus_context)
            file.write(f"{cl} - {llm_term}\n")


def string_if_not_empty(param: list) -> Union[None, str]:
    """
    Return a string if the list is not empty
    :param param: List
    :return: None if the list is empty, the string otherwise
    """
    if param == "None":
        param = []
    if param and len(param) > 0:
        l = [x for x in param if isinstance(x, float) and ~np.isnan(x) or not isinstance(x, float) and x != None]
        return "; ".join(l)
    return "no available"


@click.command(
    "ea-database", short_help="Create a database from big expression atlas files"
)
@click.option("--ea-folder", help="Expression Atlas folder", required=True)
@click.option("--ea-cl-catalog", help="Expression Atlas cell line catalog", required=True)
@click.option(
    "--output",
    help="Output file with the database",
    required=True,
    type=click.Path(exists=False),
    default="ea_database.db",
)
def ea_create_database(ea_folder: str, ea_cl_catalog: str, output: str) -> None:
    """
    The following function creates a database of celllines file from expression atlas experiments with the following information:
    each cell line will contain the following information:
    - cell line name
    - organism
    - organism part
    - age
    - developmental stage
    - sex
    - ancestry category
    - disease

    :param ea_folder: Expression Atlas folder
    :param ea_cl_catalog: Expression Atlas cell line catalog, this is a list of cell lines curated by expression atlas.
    :param output: Output file with the database
    :return:
    """

    ea_files = glob.glob(ea_folder + "/**/*.tsv", recursive=True)

    cell_lines_dict = {}
    for file in ea_files:
        # read tab-delimited file
        data = pd.read_csv(file, sep="\t")

        # remove duplicates
        data = data.drop_duplicates(
            subset=[
                "Sample Characteristic[organism]",
                "Sample Characteristic[organism part]",
                "Sample Characteristic[cell line]",
                "Sample Characteristic[disease]",
            ]
        )
        columns_data = list(data.columns)

        # add to dictionary with cell line as key
        for i, row in data.iterrows():
            cell_line = row["Sample Characteristic[cell line]"]
            if cell_line not in cell_lines_dict:
                cell_lines_dict[cell_line] = {}

                cell_lines_dict[cell_line]["organism"] = []
                cell_lines_dict[cell_line]["organism"].append(
                    row["Sample Characteristic[organism]"]
                )
                cell_lines_dict[cell_line]["organism part"] = []
                cell_lines_dict[cell_line]["organism part"].append(
                    row["Sample Characteristic[organism part]"]
                )
                cell_lines_dict[cell_line]["disease"] = []
                cell_lines_dict[cell_line]["disease"].append(
                    row["Sample Characteristic[disease]"]
                )

                # check if the other fields are present
                cell_lines_dict[cell_line]["age"] = []
                if "Sample Characteristic[age]" in columns_data:
                    cell_lines_dict[cell_line]["age"].append(
                        row["Sample Characteristic[age]"]
                    )

                cell_lines_dict[cell_line]["developmental stage"] = []
                if "Sample Characteristic[developmental stage]" in columns_data:
                    cell_lines_dict[cell_line]["developmental stage"].append(
                        row["Sample Characteristic[developmental stage]"]
                    )

                cell_lines_dict[cell_line]["sex"] = []
                if "Sample Characteristic[sex]" in columns_data:
                    cell_lines_dict[cell_line]["sex"].append(
                        row["Sample Characteristic[sex]"]
                    )

                cell_lines_dict[cell_line]["ancestry category"] = []
                if "Sample Characteristic[ancestry category]" in columns_data:
                    cell_lines_dict[cell_line]["ancestry category"].append(
                        row["Sample Characteristic[ancestry category]"]
                    )
            else:
                # check that all the fields are the same, if not raise error:
                if (
                    cell_lines_dict[cell_line]["organism"]
                    != row["Sample Characteristic[organism]"]
                ):
                    print(f"Organism is different for cell line {cell_line}")
                if (
                    cell_lines_dict[cell_line]["organism part"]
                    != row["Sample Characteristic[organism part]"]
                ):
                    print(f"Organism part is different for cell line {cell_line}")
                if (
                    row["Sample Characteristic[disease]"]
                    not in cell_lines_dict[cell_line]["disease"]
                ):
                    cell_lines_dict[cell_line]["disease"].append(
                        row["Sample Characteristic[disease]"]
                    )
                    print(
                        f"Disease is different for cell line {cell_line} - values are {cell_lines_dict[cell_line]['disease']} and {row['Sample Characteristic[disease]']}"
                    )

                if (
                    "Sample Characteristic[age]" in columns_data
                    and row["Sample Characteristic[age]"]
                    not in cell_lines_dict[cell_line]["age"]
                ):
                    cell_lines_dict[cell_line]["age"].append(
                        row["Sample Characteristic[age]"]
                    )
                    print(f"Age is different for cell line {cell_line}")

                if (
                    "Sample Characteristic[developmental stage]" in columns_data
                    and row["Sample Characteristic[developmental stage]"]
                    not in cell_lines_dict[cell_line]["developmental stage"]
                ):
                    cell_lines_dict[cell_line]["developmental stage"].append(
                        row["Sample Characteristic[developmental stage]"]
                    )
                    print(f"Developmental stage is different for cell line {cell_line}")

                if (
                    "Sample Characteristic[sex]" in columns_data
                    and row["Sample Characteristic[sex]"]
                    not in cell_lines_dict[cell_line]["sex"]
                ):
                    cell_lines_dict[cell_line]["sex"].append(
                        row["Sample Characteristic[sex]"]
                    )

                if (
                    "Sample Characteristic[ancestry category]" in columns_data
                    and row["Sample Characteristic[ancestry category]"]
                    not in cell_lines_dict[cell_line]["ancestry category"]
                ):
                    cell_lines_dict[cell_line]["ancestry category"].append(
                        row["Sample Characteristic[ancestry category]"]
                    )

                print(f"Cell line {cell_line} already in database")

    # read the cell line catalog
    ae_cl_catalog = pd.read_csv(ea_cl_catalog, sep=",", header=0)

    # check if the cell lines in the catalog are in the database
    for i, row in ae_cl_catalog.iterrows():
        if row["cell line"] in cell_lines_dict:
            print(f"Cell line {row['cell line']} found in the database")
            if row["organism"] not in cell_lines_dict[row["cell line"]]["organism"]:
                cell_lines_dict[row["cell line"]]["organism"].append(row["organism"])
            if row["organism part"] not in cell_lines_dict[row["cell line"]]["organism part"]:
                cell_lines_dict[row["cell line"]]["organism part"].append(row["organism part"])
            if row["disease"] not in cell_lines_dict[row["cell line"]]["disease"]:
                cell_lines_dict[row["cell line"]]["disease"].append(row["disease"])
            if row["age"] not in cell_lines_dict[row["cell line"]]["age"]:
                cell_lines_dict[row["cell line"]]["age"].append(row["age"])
            if row["developmental stage"] not in cell_lines_dict[row["cell line"]]["developmental stage"]:
                cell_lines_dict[row["cell line"]]["developmental stage"].append(row["developmental stage"])
            if row["sex"] not in cell_lines_dict[row["cell line"]]["sex"]:
                cell_lines_dict[row["cell line"]]["sex"].append(row["sex"])
            cell_lines_dict[row["cell line"]]["synonyms"] = [row["synonyms"]]
        else:
            print(f"Cell line {row['cell line']} not found in the database")
            cell_lines_dict[row["cell line"]] = {}
            cell_lines_dict[row["cell line"]]["organism"] = [row["organism"]]
            cell_lines_dict[row["cell line"]]["organism part"] = [row["organism part"]]
            cell_lines_dict[row["cell line"]]["disease"] = [row["disease"]]
            cell_lines_dict[row["cell line"]]["age"] = [row["age"]]
            cell_lines_dict[row["cell line"]]["developmental stage"] = [row["developmental stage"]]
            cell_lines_dict[row["cell line"]]["sex"] = [row["sex"]]
            cell_lines_dict[row["cell line"]]["synonyms"] = [row["synonyms"]]


    # write the ea atlas database to file as a comma separated file.
    with open(output, "w", newline='') as file:
        # Define the CSV headers
        headers = [
            'cell line', 'organism', 'organism part',
            'disease', 'age', 'developmental stage',
            'sex', 'ancestry category', 'synonyms'
        ]

        # Write the header row
        file.write(','.join(headers) + '\n')

        for cell_line, data in cell_lines_dict.items():
            # Construct the row
            row = [
                cell_line,
                string_if_not_empty(data.get('organism')),
                string_if_not_empty(data.get('organism part')),
                string_if_not_empty(data.get('disease')),
                string_if_not_empty(data.get('age')),
                string_if_not_empty(data.get('developmental stage')),
                string_if_not_empty(data.get('sex')),
                string_if_not_empty(data.get('ancestry category', [])),
                string_if_not_empty(data.get('synonyms', []))
            ]
            # Write the row
            file.write(','.join(row) + '\n')




@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass


cli.add_command(cl_database)
cli.add_command(nlp_recommendation)
cli.add_command(ea_create_database)

if __name__ == "__main__":
    cli()
