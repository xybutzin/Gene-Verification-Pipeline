#%%
from pathlib import Path
from Bio import SeqIO
import subprocess
import sys
project_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(project_dir))

from helper.helper import get_species_strain_name, search_sra_db

fna_dir = project_dir / 'results' / 'tblastn'
res_dir = project_dir / 'results' / 'SRA_search'

#%%
# # Uncomment this block to run single file
# # fna_file = fna_dir / 'NZ_CP011487.1_1453.fna'
# fna_file = fna_dir / 'NC_000962.3_7.fna'

# species_strain = get_species_strain_name(fna_file)
# genome_accession = species_strain.get('genome')
# species = species_strain.get('species')
# strain = species_strain.get('strain')
# genome_sra_file = res_dir / genome_accession / f'{genome_accession}_rnaseq.csv'
# print(species_strain)
# print(species)
# print(strain)
# search_sra_db(species_strain, genome_sra_file)


#%%
# Uncomment this block to run all the files
fna_files = [file for file in fna_dir.glob('*.fna') if '_upstream' not in file.name]
for fna_file in fna_files:

    # get gene_id
    gene_id = fna_file.stem
    gene_dir = res_dir / 'trg_genes' / f'{gene_id}'
    gene_dir.mkdir(parents=True, exist_ok=True)

    print(f'fna_file: {fna_file.name}')
    print(f'gene_id: {gene_id}')

    species_strain = get_species_strain_name(fna_file)
    if species_strain:
        genome_accession = species_strain.get('genome')
        genome_sra_file = res_dir / genome_accession / f'{genome_accession}_rnaseq.csv'
        
        species = species_strain.get('species')
        strain = species_strain.get('strain')
    
        print(f'genome_accession: {genome_accession}')
        print(f'species: {species}')
        print(f'strain: {strain}')
    
    if not genome_sra_file.exists():
        search_sra_db(species_strain, genome_sra_file)
    else:
        print(f'genome_sra_file {genome_sra_file.name}. Skip searching.')
    
    print('*' * 40)
     

    


    
    