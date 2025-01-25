import sys
from pathlib import Path

project_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(project_dir))

from helper.helper import combine_fasta_file, split_fasta_file


###################################
# # to combine single fasta to one fasta
# fasta_file_dir = project_dir / 'data' / 'raw'
# combine_fasta_file(pattern='NC*.faa', fasta_dir=fasta_file_dir, output='test_combined.faa')

##################################
# to split multiple fasta to single fasta
combined_fasta_file_dir = project_dir / 'data' / 'raw' / 'combined'
split_fasta_file(filename='test_combined.faa', fasta_dir=combined_fasta_file_dir)
