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

header_time_rep = []
header_size_rep = []
cdr3_rep = []
vdj_rep = []
for record in SeqIO.parse(input_path, 'fasta'):
    #print(record.id)
    header_time_rep.append(record.id.split('_')[1])
    header_size_rep.append(record.id.split('_')[2])
    cdr3_rep.append(record.id.split('_')[3])
    vdj_rep.append(record.id.split('_')[4])
    
seq_v_rep = [cdr3_rep[i] + '_' + vdj_rep[i] for i in range(len(cdr3_rep))]

lines_rep = zip(header_time_rep, seq_v_rep)
        
with open('seq_rep.csv','w') as out:
    out.write("{},{},{}\n".format( str('#'), str('tp'), str('seq_v_rep')))
    for i,j in enumerate(lines_rep):
        out.write("{},{},{}\n".format( i, j[0], j[1]))