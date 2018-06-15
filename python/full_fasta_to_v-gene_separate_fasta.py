#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 21:01:10 2018

@author: jordanzehr
"""


# In[1]:


from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(
    description='Extracts the different V genes from an unaligned fasta file')
parser.add_argument(
    '-o', '--output',
    metavar='OUTPUT',
    type=str,
    help='Directory to output results',
    required=True
)
parser.add_argument(
    '-f', '--file',
    metavar='FILE',
    type=str,
    help='Name of JSON file to parse',
    required=True
)

args = parser.parse_args()
input_path = args.file
#print(input_path)
output_directory = args.output
#print(output_directory)

# In[71]:


header = []
sequence = []
header1 = []
seq1 = []
header2 = []
seq2 = []
header3 = []
seq3 = []
header4 = []
seq4 = []
header5 = []
seq5 = []
header6 = []
seq6 = []
header7 = []
seq7 = []
headerV = []
seqV = []
bad = []


# In[72]:


def fasta_reader(title,fasta_self):
    for i,record in enumerate(SeqIO.parse(fasta_self,'fasta')):
        header.append(record.id)
        sequence.append(str(record.seq))
        if 'V1' in record.id:
            header1.append(record.id)
            seq1.append(sequence[i])
        elif 'V2' in record.id:
            header2.append(record.id)
            seq2.append(sequence[i])
        elif 'V3' in record.id:
            header3.append(record.id)
            seq3.append(sequence[i])
        elif 'V4' in record.id:
            header4.append(record.id)
            seq4.append(sequence[i])
        elif 'V5' in record.id:
            header5.append(record.id)
            seq5.append(sequence[i])
        elif 'V6' in record.id:
            header6.append(record.id)
            seq6.append(sequence[i])
        elif 'V7' in record.id:
            header7.append(record.id)
            seq7.append(sequence[i])
        elif 'V,' in record.id:
            headerV.append(record.id)
            seqV.append(sequence[i])
        else:
            bad.append(record.id)


# In[73]:


t = str(input_path.split('/')[-1] )
fasta_reader(t,input_path)


# In[42]:


def v_writer(L1,L2):
    if not L1:
        print('this list is empty')
        return
    else:
        V = str.upper(L1[0].split('_')[4].split(',')[0].split('-')[0])
        print('L1 is',V)
        print(V)
    lines = zip(L1,L2)
    #testString = str(v + ".fasta")
    #print(testString)
    out_path = str(output_directory + '/' + V + '_un.fasta')

    with open(out_path,'a') as out:
        for i, j in enumerate(lines):
            out.write('{}\n{}\n'.format('>' + j[0], j[1]))
    print(out_path)

v_writer(header1,seq1)
v_writer(header2,seq2)
v_writer(header3,seq3)
v_writer(header4,seq4)
v_writer(header5,seq5)
v_writer(header6,seq6)
v_writer(header7,seq7)
v_writer(headerV,seqV)

# In[43]:



