import json


def json_for_dashboard(input_fasta, input_json, tree, germline, output, wildcards):
    with open(input_json, 'r') as input_file:
         temp = json.load(input_file)
    CDR3_profile_coords = ['CDR3']
    FR3_profile_coords = ['FR3']

    with open(input_fasta) as file:
        fasta = file.read()
    with open(tree) as file:
        newick = file.read()

    output_dict = {
        'fasta': fasta,
        'newick': newick,
        'CDR3': CDR3_profile_coords,
        'FR3': FR3_profile_coords
    }

    with open(output, 'w') as file:
        json.dump(output_dict, file)

