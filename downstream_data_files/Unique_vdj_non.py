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

header_time_non = []
header_size_non = []
cdr3_non = []
vdj_non = []
for record in SeqIO.parse(input_path, 'fasta'):
    #print(record.id)
    header_time_non.append(record.id.split('_')[1])
    header_size_non.append(record.id.split('_')[2])
    cdr3_non.append(record.id.split('_')[3])
    vdj_non.append(record.id.split('_')[4])


seq_v_non = [cdr3_non[i] + '_' + vdj_non[i] for i in range(len(cdr3_non))]

lines_non = zip(header_time_non, seq_v_non)


with open('seq_non.csv','w') as out:
    out.write("{},{},{}\n".format( str('#'), str('tp'), str('seq_v_non')))
    for i,j in enumerate(lines_non):
        out.write("{},{},{}\n".format( i, j[0], j[1]))
        
