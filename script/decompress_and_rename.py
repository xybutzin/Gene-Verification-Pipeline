import gzip
import shutil
import os
from pathlib import Path
from helper.helper import decompress, rename

project_dir = Path(__file__).resolve().parent.parent
genome_dir = project_dir/'data'/'genome'
gz_files = genome_dir.glob('*.gz')

# Decompress the gz files
for gz_file in gz_files:
    print(f'Processing: {gz_file.name}')
    decompress(gz_file)

# Rename the files
files = [file for file in genome_dir.glob('*') if file.suffix in ['.fna', '.gff']]
for file in files:
    print(f'Renaming: {file.name}')
    # Extract the base name and modify it
    base_name = file.stem # Removes .gz suffix
    # print(f'base name: {base_name}')
    ext = file.suffix  # .gff or .fna
    prefix = '_'.join(base_name.split('_')[:2])
    new_name = f'{prefix}{ext}'
    # print(f'new name: {new_name}')
    rename(file, new_name)

