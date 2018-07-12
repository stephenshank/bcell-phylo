import json
import argparse


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

with open('data/out/%s/V%s.fasta' % (patient_id, v_gene) ) as file:
    fasta = file.read()
with open('data/out/%s/V%s.new' % (patient_id, v_gene) ) as file:
    newick = file.read()
output_dict = { 'fasta': fasta, 'newick': newick }
with open('data/out/%s/V%s.json' % (patient_id, v_gene), 'w') as file:
    json.dump(output_dict, file)
