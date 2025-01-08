from pathlib import Path
import pandas as pd
from helper.helper import download_genome

project_dir = Path(__file__).resolve().parent.parent
data_dir = project_dir/'data'/'raw'
genome_file = data_dir/'genome_accession.csv'
output_dir = project_dir/'data'/'genome'
output_dir.mkdir(parents=True, exist_ok=True)

genome_df = pd.read_csv(genome_file)
# genome_list = genome_df['genome_accession'].to_list()
accession_list = genome_df['accession_version'].to_list()

# Download the genome sequences
for i, accession in enumerate(accession_list):
    print(f'Downloading {i+1} out of {len(accession_list)} genomes: {accession}')
    download_genome(accession, output_dir)
print('Download completed')
