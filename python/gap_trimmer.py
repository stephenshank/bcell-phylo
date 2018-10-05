#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 11:44:01 2018

@author: jordanzehr
"""

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser(
    description='This code will take out the gapped regions in the codon alignment.'
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
    help='Which V-gene to build (integer from 1 to 6).',
    required=True
)

args = parser.parse_args()
patient_id = args.patient
v_gene = args.gene


infile = 'data/out/%s/V%s_codon.fasta' % (patient_id, v_gene)
codon = list(SeqIO.parse(infile, 'fasta'))



id_codon = []
seq_codon = []
for i in codon:
    id_codon.append(str(i.id))
    seq_codon.append(list(i.seq))


temp = np.array(seq_codon)

gap_fraction = np.sum(temp=='-', axis=0)/temp.shape[0]

good_sites = gap_fraction < .4

temp[:, good_sites]

list(temp[:, good_sites])

j = [''.join(row) for row in list(temp[:, good_sites])]

lines = zip(id_codon,j)
with open('data/out/%s/V%s_ungapped.fasta' % (patient_id, v_gene), 'w' ) as file:
    for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))