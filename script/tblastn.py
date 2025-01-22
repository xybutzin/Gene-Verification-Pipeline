import re
import os
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import subprocess
import sys
project_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(project_dir))

from helper.helper import extract_genome_accession_from_fasta, fetch_exact_accession, run_tblastn

faa_file_dir = project_dir / 'data' / 'raw'
faa_files = faa_file_dir.glob('*.faa')
tblastn_result_dir = project_dir / 'results' / 'tblastn'
tblastn_result_dir.mkdir(parents=True, exist_ok=True)

for faa_file in faa_files:
    gene_id = faa_file.stem
    genome_accession = extract_genome_accession_from_fasta(faa_file)
    genome_accession = fetch_exact_accession(genome_accession)
    genome_db_dir = project_dir / 'data' / 'genome' / genome_accession / genome_accession
    print(f'Running tblastn on {gene_id}')
    run_tblastn(query_faa_file=faa_file, genome_db_dir=genome_db_dir, res_dir=tblastn_result_dir)
