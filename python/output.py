import json

from Bio import SeqIO
import numpy as np


def json_for_dashboard(profile_alignment, window_trim_fas, window_trim_json, tree, germline, output, wildcards):
    window_trimmed = list(SeqIO.parse(window_trim_fas, 'fasta'))
    seq = list(SeqIO.parse(germline, 'fasta'))

    with open(window_trim_json, 'r') as input_file:
         temp = json.load(input_file)
    CDR3_profile_coords = ['CDR3']
    FR3_profile_coords = ['FR3']

    with open(window_trim_fas) as file:
        fasta = file.read()
    with open(tree) as file:
        newick = file.read()

    assert True

    output_dict = {
        'fasta': fasta,
        'newick': newick,
        'CDR3': CDR3_profile_coords,
        'FR3': FR3_profile_coords
    }

    with open(output, 'w') as file:
        json.dump(output_dict, file)

