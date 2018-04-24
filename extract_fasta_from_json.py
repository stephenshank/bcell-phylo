#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 09:49:52 2018

@author: jordanzehr
"""

import json
import itertools as it
import argparse

## change the number (3 spots) at the front of the 3 commands based on the data input

parser = argparse.ArgumentParser(
    description='Extract information from B-cell jsons.'
)
parser.add_argument(
    '-f', '--file',
    metavar='FILE',
    type=str,
    help='Name of JSON file to parse',
    required=True
)
parser.add_argument(
    '-s', '--size',
    metavar='SIZE',
    type=int,
    help='Discard entries whose size are below this threshold.',
    default=20
)

args = parser.parse_args()
input_path = args.file
size = args.size

stem, extension = input_path.split('.')
output_path = stem + '_size-' + str(size) + '_unaligned.fasta'
time = int(input_path.split('_')[0].split('/')[1])

with open(input_path, 'r') as input_file:
    data = json.load(input_file)

with open(output_path, 'w') as output_file:
    for i, item in enumerate(it.chain.from_iterable(data)):
        if int(item["size"]) > size:
            output_file.write(''.join(
                ['>seq' + str(i) + '_time-' + str(time) + '_size-' + str(item["size"]) + '_' + str(item["centroid"])[1:].replace(':','.')]
            ))
