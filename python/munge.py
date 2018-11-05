import itertools as it
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
    directory =  "data/out/" + wildcards.patient_id + "/"
    gene_fn =  directory + wildcards.v_gene + "_unique.txt"
    rearrangements_prefixes = open(gene_fn).readlines()
    rearrangement_fns = [directory + prefix.strip() +  "_unaligned_corrected_AA.fasta" for prefix in rearrangements_prefixes]
    return rearrangement_fns


def get_unique_vs(patients, clones):
    vs = []
    for patient_id in patients:
        for clone in clones:
            json_filename = 'data/input/%s_%s_clone.json' % (patient_id, clone) 
            with open(json_filename) as json_file:
                data = json.load(json_file)
            all_entries = it.chain.from_iterable(data)
            vs += [re.split(',|\*|\|', entry['tag'])[0] for entry in all_entries]
    unique_vs = list(set(vs))
    unique_vs.sort()
    with open('data/unique_vs.json', 'w') as output_file:
        json.dump(unique_vs, output_file)


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

def window_gap_trimmer(input_profile, output_fasta, output_json, wildcards):
    profile = list(SeqIO.parse(input_profile, 'fasta'))
    id_codon = [str(i.id) for i in profile]
    seq_codon = [list(str(i.seq)) for i in profile]

    CDR3_dict = {'V1':(557,564), 'V1-2':(557,564), 'V1-3':(451,458), 'V1-8':(489,496), 'V1-18':(476,483), 'V1-24':(498,505), 'V1-38-4':(679,686), 'V1-45':(432,439), 'V1-46':(583,590), 'V1-58':(581,588), 'V1-68':(6744,6751), 'V1-69':(664,671), 'V1-69-2':(681,688), 'V1-69D':(682,689), 'V2':(505,514), 'V2-5':(505,514), 'V2-10':(502,511), 'V2-26':(455,464), 'V2-70':(435,444), 'V3':(632,639), 'V3-7':(632,639), 'V3-9':(568,577), 'V3-11':(490,497), 'V3-13':(447,454), 'V3-15':(456,463), 'V3-16':(476,483), 'V3-19':(584,591), 'V3-20':(458,465), 'V3-21':(11917,11924), 'V3-22':(539,546), 'V3-23':(458,465), 'V3-25':(524,531), 'V3-29':(11917,11924), 'V3-30':(2228,2235), 'V3-30-3':(2134,2141), 'V3-32':(23127,23136), 'V3-33':(11917,11924), 'V3-35':(586,593), 'V3-38':(451,460), 'V3-43':(504,617), 'V3-43D':(699,708), 'V3-47':(286,291), 'V3-48':(622,629), 'V3-49':(678,685), 'V3-52':(655,662), 'V3-53':(481,488), 'V3-54':(585,592), 'V3-62':(6744,6751), 'V3-63':(458,467), 'V3-64':(529,536), 'V3-66':(445,452), 'V3-69':(454,461), 'V3-71':(6744,6751), 'V3-72':(541,548), 'V3-73':(978,985), 'V3-74':(471,478), 'V4':(580,587), 'V4-4':(580,587), 'V4-28':(578,585), 'V4-30-2':(292,299), 'V4-30-4':(421,438), 'V4-31':(318,325), 'V4-34':(11917,11924), 'V4-38-2':(289,294), 'V4-39':(11917,11924), 'V4-55':(658,665), 'V4-59':(2332,2339), 'V4-61':(581,588), 'V5':(301,306), 'V5-10-1':(301,306), 'V5-51':(596,603), 'V6':(777,784), 'V6-1':(777,784)}
    FR3_dict = {'V1':(443,556),'V1-2':(443,556), 'V1-3':(337,450), 'V1-8':(375,488), 'V1-18':(362,475), 'V1-24':(384,497), 'V1-38-4':(565,678), 'V1-45':(318,431), 'V1-46':(469,582), 'V1-58':(467,580), 'V1-68':(6630,6743), 'V1-69':(550,663), 'V1-69-2':(567,680), 'V1-69D':(568,681), 'V2':(391,504), 'V2-5':(391,504), 'V2-10':(388,501), 'V2-26':(341,454), 'V2-70':(321,434), 'V3':(518,631), 'V3-7':(518,631), 'V3-9':(454,567), 'V3-11':(376,489), 'V3-13':(333,446), 'V3-15':(342,455), 'V3-16':(362,475), 'V3-19':(470,583), 'V3-20':(344,457), 'V3-21':(11803,11916), 'V3-22':(425,538), 'V3-23':(344,457), 'V3-25':(410,523), 'V3-29':(11803,11916), 'V3-30':(2114,2227), 'V3-30-3':(2142,2255), 'V3-32':(23013,23126), 'V3-33':(11803,11916), 'V3-35':(472,585), 'V3-38':(337,450), 'V3-43':(618,627), 'V3-43D':(585,698), 'V3-47':(172,285), 'V3-48':(508,621), 'V3-49':(564,677), 'V3-52':(541,654), 'V3-53':(367,480), 'V3-54':(471,584), 'V3-62':(6630,6743), 'V3-63':(344,457), 'V3-64':(415,528), 'V3-66':(331,444), 'V3-69':(340,453), 'V3-71':(6630,6743), 'V3-72':(427,540), 'V3-73':(864,977), 'V3-74':(357,470), 'V4':(466,579), 'V4-4':(466,579), 'V4-28':(464,577), 'V4-30-2':(178,291), 'V4-30-4':(317,430), 'V4-31':(204,317), 'V4-34':(11803,11916), 'V4-38-2':(175,288), 'V4-39':(11803,11916), 'V4-55':(544,657), 'V4-59':(2218,2331), 'V4-61':(467,580), 'V5':(187,300), 'V5-10-1':(187,300), 'V5-51':(482,595),'V6':(663,776), 'V6-1':(777,784)}
    
    profiled_seq = None
    for seq_record in profile:
          if 'Germline' in seq_record.description:
              profiled_seq = seq_record

    profiled_np = np.array(list(str(profiled_seq.seq)), dtype='<U1')
    is_gap = profiled_np == '-'
    profile_indices = np.arange(len(is_gap))
    index_map = profile_indices[~is_gap]

    example = np.array(seq_codon, dtype='<U1')
    nongap_counts = np.sum(example != '-', axis=0)
    window = np.stack([nongap_counts[:-2], nongap_counts[1:-1], nongap_counts[2:]])
    payload = np.min(window, axis=0)
    start = np.argmax(payload >= 2)
    stop = len(payload) - np.argmax(payload[:-1] >= 2)
    reduced = example[:, start:stop+1]



    CDR3_coords = CDR3_dict[str('V'+wildcards.v_gene)]
    CDR3_profile_coords = (int(index_map[CDR3_coords[0]]- start -1), int(index_map[CDR3_coords[1]]-start-1))

    FR3_coords = FR3_dict['V'+str(wildcards.v_gene)]
    FR3_profile_coords = (int(index_map[FR3_coords[0]]-start-1), int(index_map[FR3_coords[1]]-start-1))

    j = [''.join(row) for row in list(reduced)]
    V = str('V'+wildcards.v_gene)

    output_dict = {
        'V': V,
        'CDR3': CDR3_profile_coords,
        'FR3': FR3_profile_coords
    }

    with open(output_json, 'w') as file:
        json.dump(output_dict, file)

    lines = zip(id_codon,j)
    with open(output_fasta, 'w' ) as file:
      for line in lines:
             file.write("{}{}\n{}\n".format( '>', line[0], line[1]))

