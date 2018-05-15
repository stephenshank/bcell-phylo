#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 20:29:22 2018

@author: jordanzehr
"""

from Bio import SeqIO
import argparse
 
parser = argparse.ArgumentParser(
    description='Extract information from B-cell fasta.'
)

parser.add_argument(
    '-f', '--file',
    metavar='FILE',
    type=str,
    help='Name of fasta file to parse',
    required=True
)

args = parser.parse_args()
input_path = args.file

header_time = []
header_size = []
cdr3 = []
vdj = []
for record in SeqIO.parse(input_path, 'fasta'):
    header_time.append(record.id.split('_')[1])
    header_size.append(record.id.split('_')[2])
    cdr3.append(record.id.split('_')[3])
    vdj.append(record.id.split('_')[4])


seq_v = [cdr3[i] + '_' + vdj[i] for i in range(len(cdr3))]

lines_non = zip(header_time, seq_v)

dataset = args.file.split('/')[1]
output_filename = 'data/%s/unique_vdj.csv' % dataset
with open(output_filename, 'w') as out:
    out.write("{},{},{}\n".format( str('#'), str('tp'), str('seq_v')))
    for i,j in enumerate(lines_non):
        out.write("{},{},{}\n".format( i, j[0], j[1]))
        
