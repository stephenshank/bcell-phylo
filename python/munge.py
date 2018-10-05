import re
import collections
import json
import itertools as it

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def get_clones(wildcards, clones):
    return ["data/out/" + wildcards.patient_id + "/clone_" + clone + "_unaligned.fasta" for clone in clones]


def get_rearrangements(wildcards):
    # read unique rearrangements file for each gene
    directory =  "data/out/" + wildcards.patient_id + "/"
    gene_fn =  directory + wildcards.v_gene + "_unique.txt"
    rearrangements_prefixes = open(gene_fn).readlines()
    rearrangement_fns = [directory + prefix.strip() +  "_unaligned_corrected_AA.fasta" for prefix in rearrangements_prefixes]
    return rearrangement_fns


def clone_json_to_unaligned_fasta(input, output, clone):
    with open(input, 'r') as input_file:
        data = json.load(input_file)
    bad_sequences = 0
    size = 30
    with open(output, 'w') as out:
        for i, item in enumerate(it.chain.from_iterable(data)):
            if int(item["size"]) > size:
                t = Seq(item["tag"]).split('|')[1] #this is nuc seq
                j = str(item["tag"]).split('|')[0] #this is v(d)jc
                c_list = str(item["centroid"]).split(':')
                c = c_list[-1]
            try:
                translated = str((t).translate()[1:-2])
                out.write(''.join(['>seq' + str(i) + '_time-' + str(clone) + '_size-' + str(item["size"]) + '_' + translated + '_' + str(j) + '_' +str(c)]))
            except:
                bad_sequences += 1
    if bad_sequences == 0:
        print('No bad sequencees.')
    else:
        print('%d sequences that would not translate.' % bad_sequences)


def separate_into_regions(input, output, v_gene):
    v_gene_regex = re.compile('._V(' + v_gene + ').')
    bad_searches = []
    bad_search_text = []
    v_sequences = []
    for f in input:
      current_clone = SeqIO.parse(f, 'fasta')
      for sequence in current_clone:
          regex_search = v_gene_regex.search(sequence.id)
          if not regex_search is None:
              v_gene = regex_search.group(1)
              v_sequences.append(sequence)
          else:
              bad_searches.append(sequence.id)
      SeqIO.write(v_sequences, output, 'fasta')
      print('%d searches were bad: %s' % ( len(bad_searches), bad_search_text))


def protein_and_corrected_dna(input, output_aa, output_nuc):
    records = list(SeqIO.parse(input, 'fasta'))
    k = []
    j = []
    q = []
    r = []
    s = []
    seqs_with_stops = []
    for seq_record in records:
        x = len(seq_record) % 3.0
        y = int(len(seq_record.seq.translate()))
        z = int(len(seq_record.seq.translate(to_stop=True)))
        # something might be going on here
        s = int(len(seq_record.id.split('_')[4].split('*')[0].replace('/', '-')))
        if x == 0 and y == z and s > 2:
            k.append(seq_record.description)
            q.append(seq_record.description)
            r.append(Seq(str(seq_record.seq)))
            j.append(Seq(str(seq_record.seq)).translate())
        else:
            seqs_with_stops.append('1')

    lines = zip(k,j)
    with open(output_aa, 'w') as file:
        for line in lines:
            file.write("{}{}\n{}\n".format( '>', line[0], line[1]))

    lines = zip(q,r)
    with open(output_nuc, 'w' ) as file:
        for line in lines:
            file.write("{}{}\n{}\n".format( '>', line[0], line[1]))


def protein_alignment_to_codon_alignment(protein_alignment, nucleotide_fasta, output):
    aligned_AA = list(SeqIO.parse(protein_alignment, 'fasta'))
    un_nuc = list(SeqIO.parse(nucleotide_fasta, 'fasta'))
    lines = zip(aligned_AA, un_nuc)

    for line in lines:
      u = []
      me = [str(line[1].seq[i:i+3]) for i in range(0, len(line[1]), 3)]
      this = iter(me)
      for i, j in enumerate(line[0]):
          if '-' in j:
              u.append('-'*3)
          else:
              u.append(next(this))
                
      with open(output, 'a') as out:
          out.write('>' + str(line[0].description) +'\n')
          out.write(''.join(u))
          out.write('\n')


def gap_trimmer(input, output):
    codon_file = list(SeqIO.parse(input, 'fasta'))
    id_codon = []
    seq_codon = []
    for i in codon_file:
        id_codon.append(str(i.id))
        seq_codon.append(list(i.seq))

    temp = np.array(seq_codon)
    gap_fraction = np.sum(temp=='-', axis=0)/temp.shape[0]
    good_sites = gap_fraction < .4
    j = [''.join(row) for row in list(temp[:, good_sites])]

    lines = zip(id_codon,j)
    with open(output, 'w' ) as file:
        for line in lines:
            file.write("{}{}\n{}\n".format( '>', line[0], line[1]))

