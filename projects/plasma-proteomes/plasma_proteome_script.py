import os
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
import re
import click
import requests

df = pd.DataFrame(columns=['pxd', 'doi', 'sdrf url', 'msstats url', 'method', 'samples', 'run', 'peptides'])
driver = webdriver.Chrome()

row_index = 1
row_index_msstats=1

def pubmed(paper_name):

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

def absolute(file_name, full_name, pxd_code):

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
    url_sdrf = ""

    #When the proper file is obtained the loop is interrupted by "while not". If the maximum iteration is reached, an error message is printed
    while not file_found and iteration < max_iterations:

        #First the folder with no number are checked. When it fails, the variable "number" increases, and the other folders are checked
        try:
            #In the first loop, no number is used (to avoid PXDXXXXX.0)
            if iteration == 0:
                url_options = [
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}/proteomicslfq/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}/msstatsconverter/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}/diannconverter/"
                ]
            #If it fails, in the next loop number is used and it increases up to 6
            else:
                url_options = [
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}.{number}/proteomicslfq/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}.{number}/msstatsconverter/",
                    f"http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/{pxd_code}.{number}/diannconverter/"
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

                            #If the URL is correct, it is used to obtain the URL where the sdrf file is storaged by replacing msstatsconverter/proteomicslfq/diannconverter and pipeline_info/
                            url_sdrf = url_msstats.replace("msstatsconverter/", "pipeline_info/")
                            url_sdrf = url_sdrf.replace("proteomicslfq/", "pipeline_info/")
                            url_sdrf = url_sdrf.replace("diannconverter/", "pipeline_info/") + pxd_code + ".sdrf.tsv"

                            driver.get(url_msstats)
                            file_found = True
                            break
                if file_found:
                    break
            iteration += 1
            
        except Exception as e:
            print("An error occurred:", str(e))

        number += 1

    if not file_found:
        print(file_name+"File not found within the maximum number of iterations.")

    #Both msstats and sdrf URL are stored in the dataframe
    df.loc[row_index, 'msstats url'] = msstats_url
    df.loc[row_index, 'sdrf url'] = url_sdrf

def method_run_sample(dataframe,file_list):
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
        elif 'comment[proteomics data acquisition method]' in dataframe.columns and any(item in str(dataframe.iloc[0]['comment[proteomics data acquisition method]']) for item in ['Data-Independent Adcquisition']):
            df.loc[row_index, 'method'] = "DIA"

        elif 'comment[proteomics data acquisition method]' in dataframe.columns and any(item in str(dataframe.iloc[0]['comment[proteomics data acquisition method]']) for item in ['Data-Dependent Adquisition']):
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


@click.command()
@click.option("-m", "--msstats", help="Folder where msstats files are located", required=True)
@click.option("-s", "--sdrf", help="Folder where sdrf files are located", required=True)
@click.option("-o", "--output", help="File where output is printed", required=True)
def main_script(msstats, sdrf, output):

#This is the main function that calls the rest if needed
#It uses the directory where the sdrf and msstast files are stored as input, and it creates an output file

    global row_index
    global row_index_msstats
    #For some reason, it needs these lines (global row_index, global row_index_msstats) to have access to those parameters

    folder_path = sdrf
    file_list = os.listdir(folder_path)
    msstats_folder_path = msstats
    msstats_file_list = os.listdir(msstats_folder_path)

    for file_name in file_list:
        #This first step is only to get the PXD code from the name of the file
        
        #First, it gets the directory of each file by combining the folder path and the file name
        #Then, the name of the file is obtained (PXDXXXXX) by removing the extension
        #Finally, gets the full_name by separating the name by . (This variable is stored in df)

        file_path = os.path.join(folder_path, file_name)
        dataframe = pd.read_csv(file_path, sep='\t')
        base_name = os.path.splitext(file_name)[0]
        full_name = base_name.split('.')[0]

        #Only some files contain "-", so it must be removed when present
        if '-' in base_name:
            first_part = base_name.split('-')[0]
        else:
            first_part = base_name.split('.')[0]
        df.loc[row_index, 'pxd'] = full_name
        PXD_code = first_part

        #The PXD is introduced in the URL and the API is accessed to obtain the doi
        try: 
            #Europepmc is searched by default
            europepmc_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="+ PXD_code +"&format=json"
            response = requests.get(europepmc_url)
            json_data = response.json()      
            doi_list = json_data['resultList']['result']
            doi = doi_list[0]['doi'] if doi_list else ''
            df.loc[row_index, 'doi'] = doi

        except:
            #When it can't be found in europepmc, the second step is searching in PRIDE     
            pride_url = "https://www.ebi.ac.uk/pride/ws/archive/projects/" + PXD_code

            if pride_url:
                #If the URL exists, it tries to find it
                response = requests.get(pride_url)
                json_data = response.json()
                if 'references' in json_data:
                    #First, it tries to find the doi
                    doi_list = json_data['references']
                    doi = doi_list[0]['doi'] if doi_list else ''
                    df.loc[row_index, 'doi'] = doi
                elif 'title' in json_data:
                    #If it fails, it tries to find the title of the paper, to call the function pubmed()
                    paper_name = json_data['title']
                    pubmed()
                #Otherwise, the doi is not found
                else:
                    doi = "NOT FOUND"
                    df.loc[row_index, 'doi'] = doi
            else:
                doi = "NOT FOUND"
                df.loc[row_index, 'doi'] = doi


        absolute(file_name, full_name, PXD_code)
        method_run_sample(dataframe, file_list)
        row_index += 1
    
    #In this step, the msstats files are searched to get the number of unique peptides
    for file_name in msstats_file_list:
        file_path = os.path.join(msstats_folder_path, file_name)
        dataframe = pd.read_csv(file_path, sep=',') 

        unique_values = dataframe['PeptideSequence'].nunique()
        df.loc[row_index_msstats, 'peptides'] = unique_values
        row_index_msstats += 1
        
    df.to_csv(output, index=False, sep=',')
    
if __name__ == '__main__':
    main_script() 