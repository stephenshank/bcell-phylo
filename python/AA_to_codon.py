#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 14:53:30 2018

@author: jordanzehr
"""

from itertools import cycle
from Bio import SeqIO
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

nuc = 'data/out/%s/V%s_unaligned.fasta' % (patient_id, v_gene)
aa = 'data/out/%s/V%s_AA.fasta' % (patient_id, v_gene)

aligned_AA = list(SeqIO.parse(aa, 'fasta'))
un_nuc = list(SeqIO.parse(nuc, 'fasta'))

lines = zip(aligned_AA, un_nuc)
u = []
for line in lines:
    me = [str(line[1].seq[i:i+3]) for i in range(0, len(line[1]), 3)]
    #print(me)
    counter = 0
    this = cycle(me)
    for i, j in enumerate(line[0]):
        #print(j)
        if '-' in j:
            u.append('-'*3)
        else:
            u.append(next(this))
            counter += i
    with open('data/out/%s/V%s_codon.fasta' % (patient_id, v_gene),'a') as out:
        out.write('>' + str(line[0].description) +'\n')
        out.write(''.join(u))
        out.write('\n')