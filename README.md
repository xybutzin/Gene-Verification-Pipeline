# Gene Verification Pipeline for Bacteria

This pipeline is designed to verify if a putative bacterial gene is an actual gene, focusing on claimed TRG genes (aka orphan genes). It integrates several bioinformatics tools to analyze genomic, transcriptomic, and structural data.

## Overview

This pipeline automates the verification process for bacterial genes through the following steps:

1. **Sequence Validation**: Uses tblastn to map the protein sequence to the genome and retrieve the nucleotide sequence, confirming the presence of start and stop codons.
2. **Annotation Search**: Checks genome annotation (GFF file) to determine if the sequence region is annotated.
3. **Promoter Search**: Extracts the upstream 100 bp of the sequence to predict the presence of a promoter region.(#TODO)
4. **Expression Analysis**: Searches NCBI SRA for RNA-seq experiments of the same species and strain, blasting the gene's nucleotide sequence against SRX data to check expression (#TODO: blast).
5. **Structural Prediction**: Uses AlphaFold to predict the protein structure and assess its stability.(#TODO)

This pipeline is designed for reserchers investigating the validity of claimed orphan genes (TRG genes) in bacteria.

## Project Directory Structure
```
├── bioinfo_short.yml
├── data
│   ├── genome
│   └── raw
│       ├── NC_000962.3_7.faa
│       ├── NC_000964.3_524.faa
│       └── genome_accession.csv
├── helper
│   ├── __init__.py
│   └── helper.py
├── results
│   ├── SRA_search
│   │   └── trg_genes
│   └── tblastn
├── script
│   ├── decompress_and_rename.py
│   ├── download_genome.py
│   ├── fetch_exact_accession.py
│   ├── makeblastdb.py
│   ├── parse_blast_result.py
│   ├── run_pipeline.py
│   ├── search_sra.py
│   └── tblastn.py
└── README.md
```

## Installation
Clone the repository to your local machine:
```
git clone https://github.com/xybutzin/Gene-Verification-Pipeline.git
cd Gene-Verification-Pipeline/script
```

To set up the environment and install the necessary depencies, follow these steps:

1. **Install Conda** (if not already installed):
    - Download and install Miniconda or Anaconda
2. **Set Channel Priority**
- Run the following commands to set the channel order and enable channel priority strict

```
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda config --set channel_priority strict
```
    
- Verify channel order by running:

```
conda config --show channels
```
This should display:
```
channels:
    - conda-forge
    - bioconda
    - defaults
```
3. **Create and Activate the Environment:
```
conda env create -f bioinfo_short.yml
conda activate bioinfo_short
```

## Usage
### Running the Pipeline
Run the entire pipeline by executing the main script:
```
python run_pipeline.py
```

## Workflow
### Input Files
- **Protein sequences**: Place .faa files (e.g. NC_000962.3_7.faa) in the data/raw directory.
- **Genome accession Data**: Place a .csv file (genome_accession.csv) in the data/raw directory.
    - **Format**:
        - Column 1: Bacterial name (e.g. Escherichia coli).
        - Column 2: Genome accession (e.g., GCF_003697165).

### Workflow Steps
#### Step 1: Fetch Exact Accession
- **Script**: fetch_exact_accession.py
- **Input**: genome_accession.csv
- **Process**:
    - For each genome accession in the CSV, use NCBI's **esearch** and **efetch** to retrieve the exact genome accession, including the version number.
    - Add the exact genome accession as a third column to the CSV file.
#### Step 2: Download Genome
- **Script**: download_genome.py
- **Input**: Updated genome_accession.csv with exact genome accession
- **Process**:
    - Use **ncbi-genome-download** to download the genome's FASTA and GFF files.
    - Save the files in the **data/genome** directory.
#### Step 3: Decompress and Rename
- **Script**: decompress_and_rename.py
- **Input**: Gzipped genome files in the data/genome directory
- **Process**:
    - Decompress the .gz files.
    - Rename files to follow the format <genome_accession>.fna for the FASTA file and <genome_accession>.gff for the annotation file.
#### Step 4: Make BLAST Database
- **Script**: makeblastdb.py
- **Input**: Decompressed genome FASTA files (<genome_accession>.fna)
- **Process**:
    - Use BLAST+ tools to create a database for each genome FASTA file.
#### Step 5: Run TBLASTN
- **Script**: tblastn.py
- **Input**: Protein .faa files and genome BLAST databases.
- **Process**:
    - Use **tblastn** to align protein sequences to the genome of the same strain.
    - Save the output as XML files in the **results/tblastn** directory.
#### Step 6: Parse BLAST results
- **Script**: parse_blast_result.py
- **Input**: XML files from tblastn
- **Process**:
    - Extract alignment details (e.g., % identity, % coverage, start and stop coordinates, strand, match sequences etc.).
    - Search the GFF file to check if the region is annotated.
    - Output:
        - A human-readable text file with alignment and annotation details.
        - A .fna file containing aligned nucleotide sequences.
        - A .fna file with the upstream 100 bp of the nucleotide sequence.
#### Step 7: Search SRA
- **Script**:search_sra.py
- **Input**: <protein>.faa files
- **Process**:
    - Extract exact genome accession for the strain from <protein>.faa file.
    - Use **esearch** and **efetch** to search NCBI SRA for RNA-seq experiments for the organism and strain.
    - Output a CSV file for each genome with SRX data and other relevant experiment details.


## Current Limitations
1. **Promoter Prediction**: The pipeline does not currently include a bacterial promoter prediction tool due to the lack of suitable command-line or Python-based solutions.
2. **Structural Prediction**: AlphaFold integration is pending license acquisition for use on HPC.
3. **Expression Analysis**: While the pipeline can find relevant SRX data, it cannot directly perform BLAST on the SRA database using command-line tools. Downloading SRX reads and building a local database is under consideration.

## Dependencies
- Python 3.6
- Required python libraries (see bioinfo_env.yml)
- BLAST+ tools
- Entrez Direct (EDirect)

## Contributing
Contributions are welcome! If you have suggestions for improving the pipeline or addressing its current limitations, feel free to open an issue or submit a pull request.