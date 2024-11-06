"""
This function takes the accessions from fasta file after fdrbench and add the decoy and entrapment to the uniprot accessions.
>sp|A0A087X1C5|CP2D7_HUMAN Description... -> sp|A0A087X1C5|CP2D7_HUMAN Description...
>ENTRAP_sp|A0A087X1C5|CP2D7_HUMAN Description ... -> ENTRAP_sp|ENTRAP_A0A087X1C5|ENTRAP_CP2D7_HUMAN Description...
>DECOY_sp|A0A087X1C5|CP2D7_HUMAN Description ... -> DECOY_sp|DECOY_A0A087X1C5|DECOY_CP2D7_HUMAN Description...
>DECOY_ENTRAP_sp|A0A087X1C5|CP2D7_HUMAN Description ... -> DECOY_ENTRAP_sp|DECOY_ENTRAP_A0A087X1C5|DECOY_ENTRAP_CP2D7_HUMAN Description...
"""

import os
import re
import sys


def process_id_description(line):
    id_arr = line.split('|')
    prefix = id_arr[0]
    if "ENTRAP_" in prefix and "DECOY_" in prefix:
        id_arr[1] = f'DECOY_ENTRAP_{id_arr[1]}'
        id_arr[2] = f'DECOY_ENTRAP_{id_arr[2]}'
    elif "ENTRAP_" in prefix:
        id_arr[1] = f'ENTRAP_{id_arr[1]}'
        id_arr[2] = f'ENTRAP_{id_arr[2]}'
    elif "DECOY_" in prefix:
        id_arr[1] = f'DECOY_{id_arr[1]}'
        id_arr[2] = f'DECOY_{id_arr[2]}'
    return '|'.join(id_arr)


with open('Homo-sapiens-uniprot-reviewed-contam-entrap-decoy-20241105.fasta', 'r') as file:
    f = open('output.fasta', 'w')
    # Read the contaminants file
    for line in file:
        line = line.strip()
        if line.startswith('>'):
            line = process_id_description(line)
            f.write(f'{line}\n')
        else:
            f.write(f'{line}\n')

f.close()