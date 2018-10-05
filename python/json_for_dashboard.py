import json
import argparse
import numpy as np
from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='Create JSON to be consumed by dashboard for final output.'
)

parser.add_argument(
    '-p', '--patient',
    metavar='PATIENT',
    type=str,
    help='Patient ID.',
    required=True
)

parser.add_argument(
    '-g', '--gene',
    metavar='GENE',
    type=str,
    help='Which V-gene to build (integer from 1 to 7).',
    required=True
)

args = parser.parse_args()
patient_id = args.patient
v_gene = args.gene

pro = 'data/out/%s/V%s_profile.fasta' % (patient_id, v_gene)
germ = 'data/input/Germline_nuc_V%s.fasta' % v_gene

profile = list(SeqIO.parse(pro, 'fasta'))
seq = list(SeqIO.parse(germ, 'fasta'))

profiled_seq = None
for seq_record in profile:
    if 'Gene_V' in seq_record.description:
        profiled_seq = seq_record
        print('yes')

profiled_np = np.array(list(str(profiled_seq.seq)), dtype='<U1')
is_gap = profiled_np == '-'
profile_indices = np.arange(len(is_gap))
index_map = profile_indices[~is_gap]

CDR3_dict = {'V1':(557,564), 'V2':(505,514), 'V3':(632,639), 'V4':(580,587), 'V5':(301,306), 'V6':(777,784)}
FR3_dict = {'V1':(443,556), 'V2':(391,504), 'V3':(518,631), 'V4':(466,579), 'V5':(187,300), 'V6':(663,776)}

CDR3_coords = CDR3_dict['V'+str(v_gene)]
CDR3_profile_coords = (int(index_map[CDR3_coords[0]-1]), int(index_map[CDR3_coords[1]-1]))
FR3_coords = FR3_dict['V'+str(v_gene)]
FR3_profile_coords = (int(index_map[FR3_coords[0]-1]), int(index_map[FR3_coords[1]-1]))

with open((pro)) as file:
    fasta = file.read()
with open('data/out/%s/V%s.new' % (patient_id, v_gene) ) as file:
    newick = file.read()

output_dict = {
    'fasta': fasta,
    'newick': newick,
    'CDR3': CDR3_profile_coords,
    'FR3': FR3_profile_coords
}

with open('data/out/%s/V%s.json' % (patient_id, v_gene), 'w') as file:
    json.dump(output_dict, file)
