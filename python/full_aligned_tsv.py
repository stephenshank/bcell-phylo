#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 10:24:27 2018

@author: jordanzehr
"""

# coding: utf-8

# In[1]:


from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
    description='Extract information from B-cell jsons.'
)
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
o_dir = args.output

# In[2]:


header = []
seq = []
for record in SeqIO.parse(input_path, 'fasta'):
    header.append(record.id)
    seq.append(str(record.seq))
    #print(record.id)


# In[3]:


length = []
size = []
aa_seq = []
for i in header:
    length.append(len(i.split('_')[3]))
    size.append(i.split('_')[2].split('-')[1])
    aa_seq.append(i.split('_')[3])
#print(length) # these are all ints for now, until str later
#print(len(size))
#print(len(aa_seq))
#print(size) # these are all strings
#print(aa_seq) these are all strings


# ##### seq13_time-7_size-40_TRGSRTYESGGSACD_V1,J3*02,CG3*17_250717

# In[4]:


time = []
vdj = []

##to get no reps ##
for i in header:
    n = int(i.split('_')[1].split('-')[1])
    #if not n % 2 == 0:
        #time.append(n)
    time.append(int(i.split('_')[1].split('-')[1]))
    vdj.append(i.split('_')[4])

v_gene_full = []
for i in vdj:
    v_gene_full.append(i.split(',')[0])

v_gene = []
for i in v_gene_full:
    if '*' in i:
        v_gene.append(i.split('*')[0])
    else:
        v_gene.append(i)

just_v = []
for i in v_gene:
    #print(i)
    just_v.append(i.split('-')[0])

# ### sanity check session

# In[ ]:


#print(len(time))
#print(len(v_full))
#print(len(v_gene_full))
#print(len(v_gene))
#print(len(seq))
#for i in vdj[0:2]:
#    print(i)
#for i in v_gene[0:20]:
#    print(i)


# ### we have: v full, v genes, time, aa seqs, size and length of aa seqs

# ### this will change the time ints from previous to strings. this is helpful downstream in R

# In[5]:


for index, value in enumerate(time):
    if value == int(7):
        time[index] = 'Visit:1a'
    elif value == int(8):
        time[index] = 'Visit:1b'
    elif value == int(9):
        time[index] = 'Visit:2a'
    elif value == int(10):
        time[index] = 'Visit:2b'
    elif value == int(11):
        time[index] = 'Visit:3a'
    elif value == int(12):
        time[index] = 'Visit:3b'
    else:
        next
        #time[index] = 'Visit:3b'


# In[8]:


#lines = zip(time, size, length, v_gene, v_gene_full, aa_seq)
lines = zip(time, size, length,just_v, v_gene, v_gene_full, aa_seq)

# In[9]:

p = str(o_dir + '/' + '77612_master.tsv')

#with open(p, 'w') as out:
#    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('#', 'Time', 'Size', 'CDR3 Length (AA)','V_gene','Full_V-gene', 'A.A. Seq'))
#    for i, j in enumerate(lines):
#        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( i+1, j[0], j[1], j[2], j[3],j[4],j[5]))

with open(p, 'w') as out:
    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('#', 'Time', 'Size', 'CDR3 Length (AA)','V','V_gene','Full_V-gene', 'A.A. Seq'))
    for i, j in enumerate(lines):
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( i+1, j[0], j[1], j[2], j[3],j[4],j[5],j[6]))


