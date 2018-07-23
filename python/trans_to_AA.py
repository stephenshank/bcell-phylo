#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:28:33 2018

@author: jordanzehr
"""

import argparse
from Bio.Seq import Seq
from Bio import SeqIO

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


records = list(SeqIO.parse(nuc, 'fasta'))
k= []
j = []
for seq_record in records:
    k.append(seq_record.description)
    j.append(Seq(str(seq_record.seq)).translate())


lines = zip(k,j)
with open('data/out/%s/V%s_unaligned_AA.fasta' % (patient_id, v_gene), 'w' ) as file:
    for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))