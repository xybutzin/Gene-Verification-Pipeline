import os
from pathlib import Path
import re
from Bio import SeqIO
import gffutils
from helper.helper import create_gff_db, extract_genome_accession_from_description, parse_blast_result, get_upstream_sequence

project_dir = Path(__file__).resolve().parent.parent
genome_file_dir = project_dir / 'data' / 'genome'
xml_file_dir = project_dir / 'results' / 'tblastn'
query_faa_dir = project_dir / 'data' / 'raw'

xml_files = xml_file_dir.glob('*.xml')
for xml_file in xml_files:
    # print(f'xml_file: {xml_file.name}')
    try:
        query_id = xml_file.stem
        genome_id = '_'.join(query_id.split('_')[:2])
        query_faa_file = query_faa_dir / f'{query_id}.faa'
        query = SeqIO.read(query_faa_file, 'fasta')
        
        genome_accession = extract_genome_accession_from_description(description = query.description)
        genome_file = list(genome_file_dir.glob(f'{genome_accession}.*.fna'))[0]
        # print(type(genome_file))
        # if not isinstance(genome_file, Path):
        #     genome_file = Path(genome_file)
        gff_file = genome_file.with_suffix('.gff')
        gff_db = genome_file.with_suffix('.db')
        output_txt_file = xml_file.with_suffix('.txt')
        upstream_seq_file = xml_file.with_name(xml_file.stem + '_upstream.fna')
        
        print(f'query_faa_file: {query_faa_file}')
        print(f'genome_file: {genome_file}')
        print(f'gff_file: {gff_file}')
        print(f'gff_db: {gff_db}')
        print(f'upstream_seq_file: {upstream_seq_file}')

        print(f"gff_db type: {type(gff_db)}")
        print(f"gff_file type: {type(gff_file)}")
        print(f"genome_file type: {type(genome_file)}")
        print(f"output_txt_file type: {type(output_txt_file)}")

        if not gff_db.exists():
            create_gff_db(gff_file)
        
        db = gffutils.FeatureDB(str(gff_db))
        start, end, strand = parse_blast_result(xml_file=xml_file, query_faa_file=query_faa_file, genome_file=genome_file, output_txt_file=output_txt_file)


        if start is None or end is None:
            print(f'Skipping searching for annotation for {xml_file} as start or end is None')
            continue
        
        print(f'start: {start}')
        print(f'end: {end}')
        print(f'strand: {strand}')

        with open(output_txt_file, 'a') as output_file:
            for feature in db.features_of_type('CDS'):
                # print(f'feature.start: {feature.start}')
                # print(type(feature.start))
                # print(f'feature.end: {feature.end}')
                if feature.start >= start - 50 and feature.end <= end +50:
                # if feature.start <= end and feature.end >= start:
                # if feature.start == start and feature.end == end:
                    output_file.write(f'Found annotation: from {feature.start} to {feature.end}\n')
                    output_file.write(f'ID: {feature.id}, Attributes: {feature.attributes}\n')
                    output_file.write('\n')
                    print('annotation found')
        print('-' * 40)

        with open(upstream_seq_file, 'w') as up_seq_file:
            try:
                upstream_seq = get_upstream_sequence(genome_file=genome_file, genome_id=genome_id, start=start, end=end, strand=strand, upstream_length=100)
                up_seq_file.write(f'>{query_id}_upstream_100bp\n')
                up_seq_file.write(upstream_seq)
            except ValueError as e:
                print(e)

    except Exception as e:
        print(f'An erro occured while parsing {xml_file}: {e}') 
    

# # Query CDS features matching your protein's coordinates
# start, end = 12345, 67890  # Replace with your coordinates
# for feature in db.features_of_type("CDS"):
#     if feature.start == start and feature.end == end:
#         print("Found match:")
#         print(f"ID: {feature.id}, Attributes: {feature.attributes}")


# gff_files = genome_file_dir.glob('*.gff')
# gff_column_mapping = {}
# for gff_file in gff_files:
#     print(gff_file)
#     start_col, end_col = determine_start_end_columns(gff_file)
#     print(start_col)
#     print(end_col)
#     gff_column_mapping[gff_file] = {"start_column": start_col, "end_column": end_col}







