import json
import itertools as it
import argparse
from Bio.Seq import Seq


parser = argparse.ArgumentParser(
    description='Extract information from B-cell jsons.'
)

parser.add_argument(
    '-p', '--patient',
    metavar='PATIENT',
    type=str,
    help='Patient ID.',
    required=True
)

parser.add_argument(
    '-c', '--clone',
    metavar='CLONE',
    type=int,
    help='Clone (number from 1 to 6 which specifies replicate/timepoint).',
    default=20
)

parser.add_argument(
    '-s', '--size',
    metavar='SIZE',
    type=int,
    help='Discard entries whose size are below this threshold.',
    default=20
)

args = parser.parse_args()
patient_id = args.patient
clone = args.clone
size = args.size

path_vars = ( patient_id, clone )
input_path = 'data/input/%s_%s_clone.json' % path_vars
output_path = 'data/out/%s/clone_%s_unaligned.fasta' % path_vars

with open(input_path, 'r') as input_file:
    data = json.load(input_file)

bad_sequences = 0
with open(output_path, 'w') as output_file:
    for i, item in enumerate(it.chain.from_iterable(data)):
        if int(item["size"]) > size:
            t = Seq(item["tag"]).split('|')[1] #this (should be) nuc seq
            j = str(item["tag"]).split('|')[0] #this (should be) vdj
            c_list = str(item["centroid"]).split(':')
            c = c_list[-1]
            try:
                translated = str((t).translate()[1:-2])
                output_file.write(''.join(
                ['>seq' + str(i) + '_time-' + str(clone) + '_size-' + str(item["size"])
                   + '_' + translated + '_' + str(j) + '_' +str(c) ]
                ))
            except:
                bad_sequences += 1

if bad_sequences == 0:
    print('No bad sequences.')
else:
    print('%d sequences that would not translate.' % bad_sequences)

