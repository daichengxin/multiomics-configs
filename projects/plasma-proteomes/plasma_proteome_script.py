import os
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
import re
import click
import requests
from selenium.webdriver.chrome.options import Options

chrome_options = Options()
chrome_options.add_argument('--headless')
driver = webdriver.Chrome(options=chrome_options)

row_index = 1
row_index_msstats=1

def pubmed(paper_name, df):

    #This function obtains the doi of a dataset when the doi isn't available in europepmc/PRIDE
    #I haven't been able to do it using API, so it has to use the webrowser
    #It starts by using the paper_name obtained from PRIDE, and the doi is storaged as "doi" in the dataframe df

    driver.get("https://pubmed.ncbi.nlm.nih.gov/")
    pubmed = driver.find_element(By.XPATH, "/html/body/div[2]/main/div[1]/div/form/div[1]/div[1]/div/span/input")
    pubmed.send_keys(paper_name)
    pubmed.send_keys(Keys.RETURN)

    #If the doi is not found, "Not found" is printed
    try:

        #When paper_name is typed and there is only one option, pubmed opens it automatically, so the doi can be obtained instantly
        try:
            doi_search = driver.find_element(By.XPATH, "/html/body/div[5]/main/header/div[1]/ul/li[2]/span/a")
            doi = doi_search.text
        except:
            #If there are more than one option, the first paper is selected
            #However, I found two XPATH for it, so it will try the first option and if it fails, the second one

            #First XPATH
            try:
                pubmed_articulo = driver.find_element(By.XPATH, "/html/body/main/div[9]/div[2]/section[2]/div[1]/div/article[1]/div[2]/div[1]/a")
                pubmed_articulo.click()
                #The same happens with the doi, there are two possible XPATH. It will try one and if it fails, the other one
                try:
                    doi_search = driver.find_element(By.XPATH, "/html/body/div[5]/main/header/div[1]/ul/li[2]/span/a")
                except:
                    doi_search = driver.find_element(By.XPATH, "/html/body/div[5]/main/header/div[1]/ul/li[3]/span/a")
                doi = doi_search.text

            #Second XPATH
            except:
                pubmed_articulo = driver.find_element(By.XPATH, "/html/body/main/div[9]/div[2]/section/div[1]/div/article[1]/div[2]/div[1]/a")
                pubmed_articulo.click()
                try:
                    doi_search = driver.find_element(By.XPATH, "/html/body/div[5]/main/header/div[1]/ul/li[2]/span/a")
                except:
                    doi_search = driver.find_element(By.XPATH, "/html/body/div[5]/main/header/div[1]/ul/li[3]/span/a")
                doi = doi_search.text
    except:
        pubmed_value = "NOT FOUND"
    df.loc[row_index, 'doi'] = doi

def absolute(pxd_code, full_name, file_name, df):

    #This function obtains the URL where msstats are located. 
    #Depending on the used method (DDA, DIA, TMT...) files are located in different folders, so all of them must be checked
    #In addition, when there are more than one file per PXD (multiple tissues or methods), they are named with a number (PXDXXXX.1, .2...) so that must be considered also
    #This function enters each possible folder and compares the name of the file with the pxd_code obtained previously
    #First, it enters folders with no number (proteomicslfq, msstatscoverter, diannconverter)
    #If it fails, a number is added, starting from 1, and it checks all folders again 
    #To avoid infinite loop, the maximum number is 6 (randomly selected)

    number = 0
    file_found = False
    max_iterations = 6
    iteration = 0
    msstats_url = ""
    sdrf_url = ""

    #When the proper file is obtained the loop is interrupted by "while not". If the maximum iteration is reached, an error message is printed
    while not file_found and iteration < max_iterations:

        #First the folder with no number are checked. When it fails, the variable "number" increases, and the other folders are checked
        try:
            #In the first loop, no number is used (to avoid PXDXXXXX.0)
            if iteration == 0:
                url_options = [
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}/proteomicslfq/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}/msstatsconverter/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}/diannconvert/"
                ]
            #If it fails, in the next loop number is used and it increases up to 6
            else:
                url_options = [
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}.{number}/proteomicslfq/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}.{number}/msstatsconverter/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}.{number}/diannconvert/"
                ]

            #In any case, the url is stored
            for url_msstats in url_options:
                driver.get(url_msstats)

                #Then the name of a file in the folder and the pxd_code are compared. 
                #If they are the same, it means that it is the correct folder

                #This part was difficult for me to write, and I used ChatGTP. I don't really know why (driver.find_elements(By.TAG_NAME, "a")) is used
                file_elements = driver.find_elements(By.TAG_NAME, "a")
                for file_element in file_elements:
                    if ".sdrf_openms_design_msstats_in.csv" in file_element.text:
                        msstats_absolute_file = file_element.text
                        msstats_absolute_split = msstats_absolute_file.split(".", 1)[0]
                        if full_name == msstats_absolute_split:
                            msstats_url = driver.current_url + msstats_absolute_file

                            #If the URL is correct, it is used to obtain the URL where the sdrf file is storaged by replacing msstatsconverter/proteomicslfq/diannconverter by pipeline_info/
                            sdrf_url = url_msstats.replace("msstatsconverter/", "pipeline_info/")
                            sdrf_url = sdrf_url.replace("proteomicslfq/", "pipeline_info/")
                            sdrf_url = sdrf_url.replace("diannconverter/", "pipeline_info/") + pxd_code + ".sdrf.tsv"

                            driver.get(url_msstats)
                            file_found = True
                            break
                if file_found:
                    break
            iteration += 1
            
        except Exception as e:
            print("An error occurred:", str(e))
            msstats_url="not found"
        number += 1

    if not file_found:
        print(file_name+"File not found within the maximum number of iterations.")
        msstats_url='not found'

    #Both msstats and sdrf URL are stored in the dataframe
    df.loc[row_index, 'msstats url'] = msstats_url
    df.loc[row_index, 'sdrf url'] = sdrf_url

def method_run_sample(dataframe,file_list, df):
    #This function obtains the used method (TMT, DDA...) and the number of samples and runs from the sdrf
    #TMT and iTRAQ are checked in "comment[label]" column.
    #DDA and DIA are checked in "comment[proteomics data acquisition method]"
    #If none of them is found, DDA is set by default
      
    for file_name in file_list:

        if 'TMT' in str(dataframe.iloc[0]['comment[label]']):
            df.loc[row_index, 'method'] = "TMT"

        elif 'iTRAQ' in str(dataframe.iloc[0]['comment[label]']):
            df.loc[row_index, 'method'] = "iTRAQ"    

        #For DDA and DIA, first it checks if the column 'comment[proteomics data acquisition method]' exists. If it does, it looks for "Data-Independent Acquisition" or "Data-Dependent Adquisition"
        elif 'comment[proteomics data acquisition method]' in dataframe.columns and any(item in str(dataframe.iloc[0]['comment[proteomics data acquisition method]']) for item in ['Data-Independent Acquisition']):
            df.loc[row_index, 'method'] = "DIA"

        elif 'comment[proteomics data acquisition method]' in dataframe.columns and any(item in str(dataframe.iloc[0]['comment[proteomics data acquisition method]']) for item in ['Data-Dependent Acquisition']):
            df.loc[row_index, 'method'] = "DDA"

        #If none of them is found, DDA is set by default
        else:
            df.loc[row_index, 'method'] = "DDA"
    
    #To obtain the sample number, the first column is checked.
    #First, it searches for one or more digits (\d+) at the end of the string ($) in each item (re.search(r'\d+$', item))
    #If it exists, it is converted to an integer
    #Then, the largest number is selected 

    ##sample = [int(re.search(r'Sample-(\d+)', item).group(1)) if re.search(r'Sample-(\d+)', item) else 0 for item in dataframe.iloc[:, 0]]
    sample = [int(re.search(r'\d+$', item).group()) if re.search(r'\d+$', item) else 0 for item in dataframe.iloc[:, 0]]
    largest_sample = max(sample)
    df.loc[row_index, 'samples'] = largest_sample

    #Finally, it does the same as with the samples to find numbers, but in the assay name column
    #However, this time it looks for unique values, and counts how many there are
    try:
        run = [int(re.search(r'(\d+)', item).group()) for item in dataframe['assay name']]
        run = dataframe['assay name'].nunique()
        df.loc[row_index, 'run'] = run
    except KeyError:
        df.loc[row_index, 'run'] = "error"     

def features(msstats_file_list, msstats_folder_path, df):
    #In this step, the msstats files are searched to get the number of features (not unique)
    for file_name in msstats_file_list:
        file_path = os.path.join(msstats_folder_path, file_name)
        dataframe_msstats = pd.read_csv(file_path, sep=',')

        #It gets the name of the dataset from file_name to find the proper row in df
        code = file_name.split('.')[0]
        #Then it finds the corresponding row in the DataFrame
        mask = df['pxd'] == code
        #Update the values in the 3rd column of the matching rows
        df.loc[mask, 'features'] = len(dataframe_msstats['PeptideSequence']) - 1

def peptides(peptide_file_list, peptide_folder_path, df):
    #In this step, the peptide files are searched to get the number of unique peptides
    for file_name in peptide_file_list:
        file_path = os.path.join(peptide_folder_path, file_name)
        dataframe_peptides = pd.read_csv(file_path, sep=',')

        #It gets the name of the dataset from file_name to find the proper row in df
        code = file_name.split('-peptides')[0]
        #Then it finds the corresponding row in the DataFrame
        mask = df['pxd'] == code
        #The number of unique peptides is obtained from the column 'PeptideCanonical'
        #Tha value is written in the corresponding row and column 'peptides'
        df.loc[mask, 'peptides'] = dataframe_peptides['PeptideCanonical'].nunique()

def proteins(protein_file_list, protein_folder_path, df):
    #In this step, the protein files are searched to get the number of unique proteins
    for file_name in protein_file_list:
        file_path = os.path.join(protein_folder_path, file_name)
        dataframe_protein = pd.read_csv(file_path, sep=',')

        #It gets the name of the dataset from file_name to find the proper row in df
        code = file_name.split('-proteins')[0]
        #Then it finds the corresponding row in the DataFrame
        mask = df['pxd'] == code
        #The number of unique proteins is obtained from the column 'ProteinName'
        #Tha value is written in the corresponding row and column 'proteins'
        df.loc[mask, 'proteins'] = dataframe_protein['ProteinName'].nunique()

@click.command()
@click.option("-m", "--msstats", help="Folder where msstats files are located", required=True)
@click.option("-s", "--sdrf", help="Folder where sdrf files are located", required=True)
@click.option("-m", "--peptide", help="Folder where peptide files are located", required=True)
@click.option("-s", "--protein", help="Folder where protein files are located", required=True)
@click.option("-o", "--output", help="File where output is printed", required=True)
def main_script(msstats, sdrf, peptide, protein, output):

    #This is the main function that calls the rest if needed
    #It uses the directory where the sdrf and msstast files are stored as input, and it creates an output file
    global row_index
    global row_index_msstats
    #For some reason, it needs these lines (global row_index, global row_index_msstats) to have access to those parameters

    df = pd.DataFrame(columns=['pxd', 'doi', 'sdrf url', 'msstats url', 'method', 'samples', 'run', 'features', 'peptides', 'proteins'])

    #These lines create the path to the folders where filers are storaged
    sdrf_folder_path = sdrf
    file_list = os.listdir(sdrf_folder_path)

    msstats_folder_path = msstats
    msstats_file_list = os.listdir(msstats_folder_path)

    peptide_folder_path = peptide
    peptide_file_list = os.listdir(peptide_folder_path)

    protein_folder_path = protein
    protein_file_list = os.listdir(protein_folder_path)

     #This first step is only to get the PXD code from the name of the file
    for file_name in file_list:
               
        #First, it gets the directory of each file by combining the folder path and the file name
        #Then, the name of the file is obtained (PXDXXXXX) by removing the extension
        #Finally, it gets the full_name by separating the name by . (This variable is stored in df)

        file_path = os.path.join(sdrf_folder_path, file_name)
        dataframe = pd.read_csv(file_path, sep='\t')
        base_name = os.path.splitext(file_name)[0]
        full_name = base_name.split('.')[0]

        #Only some files contain "-", so it must be removed when present
        if '-' in base_name:
            first_part = base_name.split('-')[0]
        else:
            first_part = base_name.split('.')[0]
        df.loc[row_index, 'pxd'] = full_name
        pxd_code = first_part

        #First, the pxd is introduced in europepmc. If it fails, it tries PRIDE. If the doi is not found in PRIDE, the title of the article is copied and pubmed() is called.
        #If everything fails, not found is printed

        #The PXD is introduced in the URL and the API is accessed to obtain the doi
        europepmc_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + pxd_code + "&format=json"
        response = requests.get(europepmc_url)
        json_data = response.json()
        result_data = json_data['resultList']['result']

        #If the doi is found in europepmc it is stored, otherwise, else
        if result_data and 'doi' in result_data[0] and result_data[0]['doi'] != '':
            doi = result_data[0]['doi']
            df.loc[row_index, 'doi'] = doi
        else:
            #Here it tries finding the doi in PRIDE     
            pride_url = "https://www.ebi.ac.uk/pride/ws/archive/projects/" + pxd_code
            response = requests.get(pride_url)
            #If the json exists, try finding the doi. Otherwise, doi not found
            try:
                json_data = response.json()
                doi_list = json_data['references']
                if doi_list and 'doi' in doi_list[0] and doi_list[0]['doi']:
                    doi = doi_list[0]['doi']
                    df.loc[row_index, 'doi'] = doi
                    print(pxd_code+'pride doi')
                #If the json exist, but the doi can not be found, try using the title of the article and call pubmed()
                elif 'title' in json_data:
                    paper_name = json_data['title']
                    print(pxd_code+'title')
                    pubmed(paper_name, df)
                #Otherwise, the doi is not found
                else:
                    doi = "not found"
                    df.loc[row_index, 'doi'] = doi
            except:
                doi = "not found"
                df.loc[row_index, 'doi'] = doi

        absolute(pxd_code, full_name, file_name, df)
        method_run_sample(dataframe, file_list, df)
        row_index += 1

    features(msstats_file_list, msstats_folder_path, df)
    peptides(peptide_file_list, peptide_folder_path, df)
    proteins(protein_file_list, protein_folder_path, df)
       
    #This line removes rows when the data of 'msstats url' is not found
    #this can probably be done in an easier way

    df=df.drop(df[df['msstats url'] == 'not found'].index)
     
    df.to_csv(output, index=False, sep=',')

if __name__ == '__main__':
    main_script()