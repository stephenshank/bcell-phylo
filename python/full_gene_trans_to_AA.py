#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:28:33 2018

@author: jordanzehr
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import collections

parser = argparse.ArgumentParser(
    description='Read in full gene fasta, filter out bad seqs with stop codons and ones that arent multiple of 3s, then separate the genes into all the rearrangements and make corrected AA and nuc unaligned fasta files.'
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
k = []
j = []
q = []
r = []
s = []
seqs_with_stops = []
for seq_record in records:
    x = len(seq_record)%3.0
    y = int(len(seq_record.seq.translate()))
    z = int(len(seq_record.seq.translate(to_stop=True)))
    s = int(len(seq_record.id.split('_')[4].split(',')[0].replace('/', '-')))
    if x == 0 and y == z and s > 2:
        k.append(seq_record.description)
        q.append(seq_record.description)
        r.append(Seq(str(seq_record.seq)))
        j.append(Seq(str(seq_record.seq)).translate())
    else:
        seqs_with_stops.append('1')


lines = zip(k,j)
with open('data/out/%s/V%s_unaligned_corrected_AA.fasta' % (patient_id, v_gene), 'w' ) as file:
    for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))
        
lines = zip(q,r)
with open('data/out/%s/V%s_unaligned_corrected_nuc.fasta' % (patient_id, v_gene), 'w' ) as file:
    for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))
        
###################################################################
## dont know if i need this for anything else ######################          
        
####### this is now going to separate into the rearrangements #####

infile = 'data/out/%s/V%s_unaligned_corrected_nuc.fasta' % (patient_id, v_gene)
records = list(SeqIO.parse(infile, "fasta"))

v = []
ids = []
nuc_seqs = []
for seq_record in records:
    if '/' in seq_record.id:
        ## this will now NOT grab different alleles
        temp = seq_record.id.split('_')[4].split(',')[0].split('*')[0].replace('/', '-')
        #temp = seq_record.id.split('_')[4].split(',')[0].replace('/', '-')
        v.append(temp)
        temp2 = seq_record.id.replace('/', '-')
        ids.append(str(temp2))
        nuc_seqs.append(str(seq_record.seq))
    else:
        v.append(seq_record.id.split('_')[4].split(',')[0].split('*')[0])
        ids.append(str(seq_record.id))
        nuc_seqs.append(str(seq_record.seq))
        
v_count = collections.Counter(v)
uni_v = []
for i, j in v_count.items():
    uni_v.append(i)




##### now just writing the correct files for our individual rearrangements ###
for e in uni_v:
    for pos, i in enumerate(ids):
        if e in i:
            outfile = 'data/out/%s/%s_unaligned_corrected_nuc.fasta' % (patient_id,e)
            with open(outfile,'a') as out:
                out.write("{}{}\n{}\n".format('>',str(i),str(nuc_seqs[pos])))
        else:
            next

##### same loop, just adding the translated Amino Acid part to make this fasta ### 
for e in uni_v:
    for pos, i in enumerate(ids):
        if e in i:
            outfile = 'data/out/%s/%s_unaligned_corrected_AA.fasta' % (patient_id,e)
            with open(outfile,'a') as out:
                temp = Seq(nuc_seqs[pos])
                trans = str(temp.translate())
                out.write("{}{}\n{}\n".format('>',str(i),str(trans)))
        else:
            next
            

outfile = 'data/out/%s/V%s_unique.txt' % (patient_id,v_gene)
with open(outfile, 'w') as out:
    for i in uni_v:
        out.write("%s\n" % i)
         
         
         
         
         
         
         
         