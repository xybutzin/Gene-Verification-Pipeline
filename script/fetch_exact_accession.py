import os
from pathlib import Path
import pandas as pd
import subprocess
import sys
project_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(project_dir))

from helper.helper import fetch_exact_accession

data_dir = project_dir/'data'/'raw'
genome_file = data_dir/'genome_accession.csv'

genome_df = pd.read_csv(genome_file)
genome_list = genome_df['genome_accession'].to_list()
# print(genome_list)

genome_accession_list = []
for genome in genome_list:
    try:
        genome_accession = fetch_exact_accession(genome)
        genome_accession_list.append(genome_accession)
        
    except Exception as e:
        print(f"Error: {e}")
        genome_accession_list.append(None)
    

genome_df['accession_version'] = genome_accession_list
genome_df.to_csv(genome_file, index=False)