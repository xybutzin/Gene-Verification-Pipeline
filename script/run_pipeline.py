import subprocess

pipeline = [
    'fetch_exact_accession.py',
    'download_genome.py',
    'decompress_and_rename.py',
    'makeblastdb.py',
    'tblastn.py',
    'parse_blast_result.py',
    'search_sra.py'
]

for script in pipeline:
    subprocess.run(['python', script], check=True)