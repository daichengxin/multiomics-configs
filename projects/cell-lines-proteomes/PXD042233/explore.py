import os
import urllib
from contextlib import closing

import pandas as pd
import requests
from tqdm import tqdm

pd.options.display.max_columns = 80

# Check if 'pride_metadata.csv' is in the current folder, if not download it
if 'pride_metadata.csv' not in os.listdir():
    response = requests.get('https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD042233/pride_metadata.csv', stream=True)
    with open("pride_metadata.csv", "wb") as handle:
        for data in tqdm(response.iter_content()):
            handle.write(data)

df = pd.read_csv('pride_metadata.csv', index_col=0)
df.head()