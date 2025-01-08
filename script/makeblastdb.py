from pathlib import Path
from helper.helper import makeblastdb

project_dir = Path(__file__).resolve().parent.parent
genome_dir = project_dir / 'data' / 'genome'

genome_fna_files = genome_dir.glob('*.fna')

for genome_fna_file in genome_fna_files:
    print(f'Make blast database for {genome_fna_file.name}')
    makeblastdb(genome_fna_file)