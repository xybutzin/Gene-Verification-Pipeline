import sys
from pathlib import Path
import re
import gzip
import shutil
import subprocess
from Bio import SeqIO
import pandas as pd
import gffutils
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from xml.etree import ElementTree as ET


def combine_fasta_file(pattern=None, fasta_dir=Path.cwd(), output='combined'):
    '''
    combine multiple fasta (.fasta, .fna, or .fas) files into a single file.
    Default input folder is the current working directory
    '''

    fasta_dir = Path(fasta_dir)
    if pattern:
        fasta_files = fasta_dir.glob(pattern)
    else:
        fasta_files = fasta_dir.glob('*')

    try: 
        first_fasta_file = next(fasta_files)
    except StopIteration:
        print("No matching files found in the directory")
        return None
    
    f_ext = first_fasta_file.suffix
    if not Path(output).suffix:
        output = fasta_dir/f'{output}{f_ext}'
    else:
        output = fasta_dir/f'{output}'
    print(output)

    with open(output, 'w') as out_file:
        # Write the first file
        with open(first_fasta_file, 'r') as in_file:
            print(f'processing: {first_fasta_file.name}')
            for line in in_file:
                out_file.write(line)
            out_file.write('\n')

        # Write the remaining files
        for fasta_file in fasta_files:
            print(f'processing: {fasta_file.name}')
            with open(fasta_file, 'r') as input:
                for line in input:
                    out_file.write(line)
                out_file.write('\n')

def split_fasta_file(filename, fasta_dir= Path.cwd()):

    if not isinstance(fasta_dir, Path):
        fasta_dir = Path(fasta_dir)
    combined_file = fasta_dir / f'{filename}'
    f_ext = combined_file.suffix
    print(f_ext)

    for record in SeqIO.parse(combined_file, 'fasta'):
        id = record.id
        description = record.description
        sequence = str(record.seq)
        single_file = fasta_dir/f'{id}{f_ext}'
        print(f'Creating {single_file.name}')
        with open(single_file, 'w') as f:
            f.write(f'>{description}\n')
            f.write(sequence)


# Function to run esearch, efetch, grep pipeline to get exact accession
def fetch_exact_accession(query_accession):
    '''
    Take query_accession (GCF_003697165) and 
    output genome accession with version (GCF_003697165.2)
    '''
    command = (
        f'esearch -db assembly -query "{query_accession}" | '
        'efetch -format docsum | '
        'grep -oP "(?<=<RefSeq>).*?(?=</RefSeq>)" | head -1'
    )

    # Execute the command pipeline
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    # Check for errors
    if process.returncode != 0:
        raise RuntimeError(f"Command failed with error: {stderr.decode('utf-8')}")
    
    # Decode and strip the output
    return stdout.decode('utf-8').strip()



# Function to run ncbi-genome-download
def download_genome(accession, output_dir):
    '''
    Input is genome accession with version (GCF_003697165.2). Downloads the genome 
    in fasta format (nt seq) and gff format (annotation) in the output_dir
    '''

    command = [
        "ncbi-genome-download",
        "bacteria",
        "-s", "refseq",
        "--assembly-accessions", accession,
        "-F", "gff,fasta",
        "-o", str(output_dir),
        "--flat-output"
    ]

    subprocess.run(command)


def decompress(gz_file):
    '''decompress .gz file'''
    
    if not isinstance(gz_file, Path):
        gz_file = Path(gz_file)
        
    output_file = gz_file.with_suffix('')  # removes .gz
    
    with gzip.open(gz_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    # Remove the orginial .gz file
    gz_file.unlink()


def rename(old_name, new_name):
    '''
    Rename old_name or directory.
    Args:
        old_name: Path to the file or dir to rename, can be a string or Path object.
        new_name: New name for the file or dir (string)
    '''
    if not isinstance(old_name, Path):
        old_name = Path(old_name)
    file_dir = old_name.resolve().parent
    new_name = file_dir/new_name
    try:
        old_name.rename(new_name)
    except Exception as e:
        print(f'Error renaming file: {e}')

    
def makeblastdb(genome_fna_file):
    '''
    Create a BLAST database for the genome
    args:
        input is the genome.fna file (i.e GCF_003697165.2.fna)
        output db is stored in the directory with the genome_accession name
    '''

    if not isinstance(genome_fna_file, Path):
        genome_fna_file = Path(genome_fna_file)

    # Create a directory to store the db
    genome_dir = genome_fna_file.resolve().parent
    base_name = genome_fna_file.stem
    db_dir = genome_dir / base_name
    db_dir.mkdir(parents=False, exist_ok=True)

    # Run makeblastdb and save output in the directory created
    try:
        subprocess.run(['makeblastdb', 
                        '-in', genome_fna_file, 
                        '-dbtype', 'nucl',
                        '-out', db_dir / base_name,
                        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running makeblastdb: {e}")
    

def extract_genome_accession_from_fasta(fasta_file):
    '''
    Extracts the genome accession (no version number, ie GCF_003697165) 
    from the head(description) of protein or nt sequence (fasta/fna/faa file)
    Args:
        input fasta_file (e.g. NC_000962.3_7.faa or NC_000962.3_7.fna)
    '''
    record = SeqIO.read(fasta_file, 'fasta')
    description = record.description
    match = re.search(r"genome:GCF_(\d+)", description)
    if match:
        return f'GCF_{match.group(1)}'
    else:
        return None

def extract_genome_accession_from_description(description):
    '''
    Extracts the genome accession (no version number, ie GCF_003697165) 
    from the head(description) of protein or nt sequence (fasta/fna/faa file)
    Args:
        input is header/description of a fasta_file (e.g. NC_000962.3_7.faa or NC_000962.3_7.fna)
    '''
    match = re.search(r"genome:GCF_(\d+)", description)
    if match:
        return f'GCF_{match.group(1)}'
    else:
        return None   

def get_accession_with_version_from_csv(genome_csv, accession):
    '''Fetches the accession version for a given GCF_number from the genome csv file'''
    genome_df = pd.read_csv(genome_csv)
    try:
        return genome_df.loc[genome_df['genome_accession'] == accession, 'accession_version'].iloc[0]
    except IndexError:
        print(f"Accession {accession} not found in {genome_file}")
        return None


def run_tblastn(query_faa_file, genome_db_dir, res_dir):
    '''
    Performs tblastn for query protein sequence in the given FASTA file.
    results stored in res_dir
    '''

    record = SeqIO.read(query_faa_file, 'fasta')
    # Extract information from the record
    sequence_id = record.id
    sequence = str(record.seq)

    if not isinstance(genome_db_dir, Path):
        genome_db_dir = Path(genome_db_dir)
    genome_accession = genome_db_dir.name
    # print(f'genome accession: {genome_accession}')
    
    if not genome_accession:
        print(f"Error getting genome accession. Skipping {sequence_id}.")
    
    result_file = res_dir / f'{sequence_id}.xml'
    # print(result_file)

    tblastn_cmd = [
        'tblastn',
        '-query', query_faa_file,
        '-db', genome_db_dir,
        '-out', result_file,
        '-outfmt', '5',
        '-evalue', '1e-5',
        '-num_threads', '4'
    ]

    try:
        # print(f'Running tblastn on {sequence_id}')
        subprocess.run(tblastn_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during tblastn with {sequence_id}: {e}")

def run_tblastn_on_all_proteins(multi_faa_file):
    '''Performs tblastn for each protein sequence in the given FASTA file containing multiple fasta seqs.'''

    pass

def create_gff_db(gff_file):

    if not isinstance(gff_file, Path):
        gff_file = Path(gff_file)

    db_file = gff_file.with_suffix('.db')
    print(db_file)

    # gffutils argments takes str, not pathlib.PosixPath
    gffutils.create_db(
        str(gff_file),
        dbfn=str(db_file),
        keep_order=True,
        merge_strategy='merge',
        verbose=False
    )
    # print('database created')

def parse_blast_result(xml_file, query_faa_file, genome_file, output_txt_file):
    '''
    parse blast result in xml file and write to a txt file
    write aligned nt seq to .fna file
    '''
    # print(f"xml_file type: {type(xml_file)}")
    # print(f"query_faa_file type: {type(query_faa_file)}")
    # print(f"genome_file type: {type(genome_file)}")
    # print(f"output_txt_file type: {type(output_txt_file)}")

    if not isinstance(query_faa_file, Path):
        query_faa_file = Path(query_faa_file)
    if not isinstance(xml_file, Path):
        xml_file = Path(xml_file)

    query_id = query_faa_file.stem    # NC_000962.3_7
    query_description = SeqIO.read(query_faa_file, 'fasta').description
    
    genome_id = '_'.join(query_id.split('_')[:2])  # NC_000962.3
    genome_sequences = {record.id: record.seq for record in SeqIO.parse(genome_file, 'fasta')}
    genome_sequence = genome_sequences[genome_id]
    # query_fna_file = query_faa_file.with_suffix('.fna')
    query_fna_file = xml_file.with_name(query_faa_file.stem + '.fna')

    try:
        with open(xml_file) as blast_result_handle, \
            open(output_txt_file, 'w') as output_file, \
            open(query_fna_file, 'w') as nt_file:
            
            blast_records = NCBIXML.parse(blast_result_handle)
            for record in blast_records:
                output_file.write(f'Query: {record.query}\n')   # for ours, there is just a single query == one record
                if record.alignments:
                    for alignment in record.alignments:
                        # print(f"\nHit: {alignment.hit_def}")
                        # print(f"Length: {alignment.length}")


                        for hsp in alignment.hsps:
                            identity = (hsp.identities / hsp.align_length) * 100
                            strand = "Forward" if hsp.frame[1] > 0 else "Reverse Complement"
                            
                            query_seq = hsp.query
                            match_seq = hsp.match
                            subject_seq = hsp.sbjct
                        
                            start = hsp.sbjct_start -1
                            end = hsp.sbjct_end

                            coverage = hsp.align_length / ((end - start) / 3 ) * 100

                            if strand == 'Forward':
                                aligned_nucleotide_seq = genome_sequence[start:end]
                            else:
                                aligned_nucleotide_seq = genome_sequence[start:end].reverse_complement()

                            # Write alignment results to the .txt file
                            output_file.write(f'Alignment: {alignment.title}\n')
                            output_file.write(f'E-value: {hsp.expect}\n')
                            output_file.write(f'% Identify: {identity:.2f}\n')
                            output_file.write(f'% Coverage: {coverage:.2f}\n')
                            output_file.write(f'Query Sequence:  {query_seq}\n')
                            output_file.write(f'Match Sequence:  {match_seq}\n')
                            output_file.write(f'Subject seuqnce: {subject_seq}\n')
                            output_file.write(f'Aligned Nucleotide Sequence: {aligned_nucleotide_seq}\n')
                            output_file.write(f'Strand: {strand}\n')
                            output_file.write(f'Start: {start + 1}, End: {end}\n')
                            output_file.write('-' * 40 + '\n')

                            # Write the aligned nt seq to a fasta (.fna) file
                            nt_file.write(f'>{query_id} start:{start+1} end:{end} strand:{strand} {query_description}\n{aligned_nucleotide_seq}')
                            return(start + 1, end, strand)
                else:
                    output_file.write("No alignments found.\n")
                    output_file.write('-' * 40 + '\n')
    except Exception as e:
        print(f"An error occured: {e}")



def determine_start_end_columns(gff_file):
    '''
    function to determine columns for the start and end of CDS features in the gff file
    '''
    if isinstance(gff_file, Path):
        gff_file = str(gff_file)

    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            if not 'CDS' in line:
                continue
            columns = line.strip().split('\t')
            # print(columns)
            index = columns.index('CDS')
            start_col = index + 1 # python 0-based index
            end_col = index + 2 # python 0-based index
            return start_col, end_col
        

def get_upstream_sequence(genome_file, genome_id, start, end, strand, upstream_length=100):
    '''
    Extracts the upstream sequence of a given length from a reference point, save as _upstream.fna
    '''

    # parse the genome file to find the sequence with the specified genome_id
    for record in SeqIO.parse(genome_file, 'fasta'):
        if record.id == genome_id:
            genome_sequence = record.seq

            if strand == 'Forward':
                reference_index = start - 1  # convert 1-based reference point to 0-based index

                # Ensure upstream sequence is withine bounds
                if reference_index < upstream_length:
                    print(f'Reference point is too close to the start of the sequence')
                    upstream_sequence = genome_sequence[:reference_index]  # take all bases before reference index
                else:
                    upstream_sequence = genome_sequence[reference_index - upstream_length:reference_index]
                

            else:
                reference_index = end

                if reference_index + upstream_length > len(genome_sequence):
                    print(f'Reference point is too close to the end of the sequence')
                    upstream_sequence = genome_sequence[reference_index:].reverse_complement()
                else: 
                    upstream_sequence = genome_sequence[reference_index:reference_index + upstream_length].reverse_complement()
            return str(upstream_sequence)
    raise ValueError(f'Genome ID {genome_id} not found in the {genome_file} file')


def get_species_strain_name(fna_file):
    '''
    Use eserach and efetch to retrieve the species and strain name
    '''
    try:
        record = SeqIO.read(fna_file, 'fasta')
    except Exception as e:
        print(f'Error when reading {fna_file.name}: {e}')
        return None

    sequence_id = record.id
    description = record.description

    genome_accession = extract_genome_accession_from_description(description)
    
    if genome_accession:
        genome_accession = fetch_exact_accession(genome_accession)
        # genome_accession = 'GCF_000195955.2'
        # print(genome_accession)

        esearch_cmd = [
        'esearch',
        '-db', 'assembly',
        '-query', genome_accession
        ]
        efetch_cmd = ['efetch', '-format', 'docsum']

        esearch_process = subprocess.Popen(esearch_cmd, stdout=subprocess.PIPE)
        efetch_process = subprocess.Popen(efetch_cmd, stdin=esearch_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        xml_output, error = efetch_process.communicate()

        # Decode the xml_output to string
        xml_output = xml_output.decode('utf-8')

        if error:
            print(f'Error: {error.decode("utf-8")}')
            return None
        
        # print(xml_output)
        species_strain = {}

        tree = ET.ElementTree(ET.fromstring(xml_output))
        root = tree.getroot()

        for doc_summary in root.findall('DocumentSummary'):
            
            organism = doc_summary.find('Organism').text
            species = doc_summary.find('SpeciesName').text
            species_strain['genome'] = genome_accession
            species_strain['organism'] = organism
            species_strain['species'] = species
            # print(species)

            biosource = doc_summary.find('Biosource')
            infraspecies_list = biosource.find('InfraspeciesList')

            for infraspecie in infraspecies_list.findall('Infraspecie'):
                sub_type = infraspecie.find('Sub_type').text
                if sub_type == 'strain':
                    strain = infraspecie.find('Sub_value').text
                    species_strain['strain'] = strain
    
        return species_strain
        
    else:
        print(f"No match found for genome_accession for {sequence_id}")
        return None
    

def search_sra_db(species_strain, csv_output):
    '''
    Takes species and strain as input query and search in the SRA database,
    outputs a list of SRR that are RNASeq for the species and strain
    '''
    genome_accession = species_strain.get('genome')
    species = species_strain.get('species')
    strain = species_strain.get('strain')

    if not isinstance(csv_output, Path):
        csv_output = Path(csv_output)
    csv_output.parent.mkdir(parents=True, exist_ok=True)
    print(f'csv_output: {csv_output}')

    esearch_cmd = [
        'esearch',
        '-db', 'sra',
        '-query', f'"{species}"[Organism] AND "{strain}"[All Fields] AND "RNA-seq"[strategy]'
    ]

    efetch_cmd = [
        'efetch',
        '-format', 'runinfo'
    ]

    with open(csv_output, 'w') as outfile:
        esearch_process = subprocess.Popen(esearch_cmd, stdout=subprocess.PIPE)
        efetch_process = subprocess.Popen(efetch_cmd, stdin=esearch_process.stdout, stdout=outfile)
        esearch_process.stdout.close()  # Allow esearch to receive a SIGPIPE if efetch exisits)
        efetch_process.communicate()
