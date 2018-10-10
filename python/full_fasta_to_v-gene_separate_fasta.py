import argparse
import re

from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Extracts the different V genes from a given patient'
)

parser.add_argument(
    '-p', '--patient',
    metavar='PATIENT',
    type=str,
    help='Patient ID.',
    required=True
)

args = parser.parse_args()
patient_id = args.patient

all_v_sequences = [[] for i in range(7)]
v_gene_regex = re.compile('._V(\d).')
bad_searches = []
for clone in [1, 2, 3, 4, 5, 6]:
    path_vars = ( patient_id, clone )
    input_path = 'data/out/%s/clone_%d_unaligned.fasta' % path_vars
    current_clone = SeqIO.parse(input_path, 'fasta')
    for sequence in current_clone:
        regex_search = v_gene_regex.search(sequence.id)
        if not regex_search is None:
            v_gene = int(regex_search.group(1))
            sequence.id = sequence.id.replace(',', '-')
            sequence.name = ''
            sequence.description = ''
            assert not ',' in sequence.id
            all_v_sequences[v_gene-1].append(sequence)
        else:
            bad_searches.append(sequence.id)

for i, v_sequences in enumerate(all_v_sequences):
    output_path = 'data/out/%s/V%d_unaligned.fasta' % (patient_id, i+1)
    SeqIO.write(v_sequences, output_path, 'fasta')

bad_search_text = '\n'.join(bad_searches)
print('%d searches were bad: %s' % ( len(bad_searches), bad_search_text) )
