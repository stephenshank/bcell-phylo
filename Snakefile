include: "rules/common.smk"

#PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]


def debug(wildcards):
  import pdb;pdb.set_trace()
## need to make the quality adjustments to the sequences

## need to add all the patient ids
#PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
PATIENT_IDS = ["28729"]
CLONES = ["1", "2", "3", "4", "5", "6"]
### this set contains OR genes which arent in tar ball. not on IMGT.
# GENES = ["1", "1-18", "1-2", "1-24", "1-3", "1-46", "1-69", "1-69-2", "1-69D", "1-8", "1-45", "1-58", "1-OR15-1", "1-OR15-9", "2", "2-5", "2-70", "2-26", "2-70D", "3", "3-11", "3-13", "3-15", "3-20", "3-21", "3-23", "3-30", "3-30-3", "3-33", "3-43", "3-48", "3-49", "3-53", "3-64", "3-66", "3-7", "3-72", "3-74", "3-9", "3-OR16-13", "3-43D", "3-73", "3-OR163-OR16-9", "4", "4-30", "4-30-2", "4-30-4", "4-31", "4-34", "4-38-2", "4-39", "4-4", "4-59", "4-61", "4-28", "4-OR15-8", "5", "5-51", "5-10-1", "6", "6-1"]
### this is just the base set 
#GENES = ["1", "2", "3", "4", "5", "6"]
### had to take out V2-70D, no 4-30, v4-28 only one sequene so took it out
GENES = ["1", "1-18", "1-2", "1-24", "1-3", "1-46", "1-69", "1-69-2", "1-69D", "1-8", "1-45", "1-58", "2", "2-5", "2-70", "2-26", "3", "3-11", "3-13", "3-15", "3-20", "3-21", "3-23", "3-30-3", "3-33", "3-43", "3-48", "3-49", "3-53", "3-64", "3-66", "3-7", "3-72", "3-9", "3-43D", "3-73", "3-74", "4", "4-30-2", "4-30-4", "4-31", "4-34", "4-38-2", "4-39", "4-4", "4-59", "4-61", "5", "5-51", "5-10-1", "6", "6-1"]
#GENES = ["1-18"]


######################################
### This input should be the absolute last 
### file to be produced from the pipeline
######################################
rule all:
  input:
  ### clone_unaligned_fasta output !!!!!! something is weird for this expand
    #expand("data/out/{patient_id}/clone_{clone}_unaligned.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ### gene_unaligned_fasta output 
    #expand("data/out/{patient_id}/V{v_gene}_unaligned.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ### unaligned_to_amino_acids output
    #expand("data/out/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ###  alignments output
    #expand("data/out/{patient_id}/V{v_gene}_aligned_AA.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ###  codon_maker output
    #expand("data/out/{patient_id}/V{v_gene}_codon.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ###  profile_alignments output
    #expand("data/out/{patient_id}/V{v_gene}_profile.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ###  gap_trimmer output
    #expand("data/out/{patient_id}/V{v_gene}_ungapped.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ###  trees output
    #expand("data/out/{patient_id}/V{v_gene}.new", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
  ###  v_gene_for_json output
    expand("data/out/{patient_id}/V{v_gene}.json", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)

######################################
### This rule reads in clone files then
### builds unaligned fasta files from them
######################################
rule clone_unaligned_fasta:
  input:
    cln="data/input/{patient_id}_{clone}_clone.json"
  ####
  output:
    cln_v="data/out/{patient_id}/clone_{clone}_unaligned.fasta"
  ####
  run:
    #import pdb;pdb.set_trace()
    import json
    import itertools as it
    from Bio.Seq import Seq

    with open(input.cln, 'r') as input_file:
      data = json.load(input_file)
    bad_sequences = 0
    size = 30
    with open(output.cln_v, 'w') as out:
      for i, item in enumerate(it.chain.from_iterable(data)):
        if int(item["size"]) > size:
          t = Seq(item["tag"]).split('|')[1] #this is nuc seq
          j = str(item["tag"]).split('|')[0] #this is v(d)jc
          c_list = str(item["centroid"]).split(':')
          c = c_list[-1]
          try:
            #print('asjhdflaksdhflaskhdflasdhflasdjhfl')
            translated = str((t).translate()[1:-2])
            #out.write('dog')
            # ]))
            out.write(''.join(['>seq' + str(i) + '_time-' + str(wildcards.clone) + '_size-' + str(item["size"]) + '_' + translated + '_' + str(j) + '_' +str(c)]))
  
          except:
            bad_sequences += 1
    if bad_sequences == 0:
      print('No bad sequencees.')
    else:
      print('%d sequences that would not translate.' % bad_sequences)

######################################
### This is where we can separate based on 'rearrangment'
### thus rule reads in the newly built fasta files
### and separates the fastas to unique V region fastas
######################################
rule gene_unaligned_fasta:
  input:
    clones=get_clones
    #grub=rules.clone_unaligned_fasta.output.cln_v
  ####
  output:
    v_gene_output="data/out/{patient_id}/V{v_gene}_unaligned.fasta"
  ####
  log:
    "logs/{patient_id}_{v_gene}_clone_to_gene.log"
  ####
  run:
    #import pdb;pdb.set_trace()
    v_gene_regex = re.compile('._V(' + wildcards.v_gene + ').')
    bad_searches = []
    bad_search_text = []
    v_sequences = []
    for f in input.clones:
      current_clone = SeqIO.parse(f, 'fasta')
      for sequence in current_clone:
          regex_search = v_gene_regex.search(sequence.id)
          if not regex_search is None:
              v_gene = regex_search.group(1)
              #v_gene = int(regex_search.group(1))
              v_sequences.append(sequence)
          else:
              bad_searches.append(sequence.id)
      SeqIO.write(v_sequences, output.v_gene_output, 'fasta')
      print('%d searches were bad: %s' % ( len(bad_searches), bad_search_text))

######################################
### we can manipulate this python script to get the individual genes ####
### This rule reads in each unaligned V region fasta file
### does a quality check on the translated seqs (seq divisible by 3, no stop codonds in seq, seq properly annotated)
### writes these new amino acid fasta files, then goes back to the original nucleotide fasta file and 
### makes a corresponding fasta file to the 'corrected' amino acid fasta
######################################
rule unaligned_to_amino_acids:
  input:
    un="data/out/{patient_id}/V{v_gene}_unaligned.fasta",
    #"python/full_gene_trans_to_AA.py"
  ####
  output:
    unaligned_aa="data/out/{patient_id}/V{v_gene}_unaligned_corrected_AA.fasta",
    un_cor_nuc="data/out/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta"
  ####
  run:
    from Bio import SeqIO
    from Bio.Seq import Seq
    import collections

    records = list(SeqIO.parse(input.un, 'fasta'))
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
    with open(output.unaligned_aa, 'w') as file:
      for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))

    lines = zip(q,r)
    with open(output.un_cor_nuc, 'w' ) as file:
      for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))
  #shell:
    #"python python/full_gene_trans_to_AA.py -p {wildcards.patient_id} -g {wildcards.v_gene}"

######################################
### This rule aligns the amino acid fasta files
######################################
rule alignments:
  input:
    rules.unaligned_to_amino_acids.output.unaligned_aa
  ####
  output:
    dog="data/out/{patient_id}/V{v_gene}_aligned_AA.fasta",
  ####
  shell:
    "mafft --amino {input} > {output}"

######################################
### This rule reads in the aligned amino acid fasta and
### the unaligned corresponding nucleotide fasta file to make a codon fasta
######################################
rule codon_maker:
  input:
    #test=debug
    aligned_aa=rules.alignments.output.dog,
    un_cor_nuc=rules.unaligned_to_amino_acids.output.un_cor_nuc
  ####
  output:
    codon="data/out/{patient_id}/V{v_gene}_codon.fasta",
  ####
  run:
    #import pdb;pdb.set_trace()
    aa_file=input.aligned_aa
    nuc_file=input.un_cor_nuc
    #import pdb;pdb.set_trace()
    aligned_AA = list(SeqIO.parse(input.aligned_aa, 'fasta'))
    un_nuc = list(SeqIO.parse(input.un_cor_nuc, 'fasta'))
    #import pdb;pdb.set_trace()  
    lines = zip(aligned_AA, un_nuc)

    for line in lines:
      u = []
      me = [str(line[1].seq[i:i+3]) for i in range(0, len(line[1]), 3)]
        #print(me)
        #counter = 0
      this = iter(me)
      for i, j in enumerate(line[0]):
            #print(j)
          if '-' in j:
              u.append('-'*3)
          else:
              u.append(next(this))
                
      with open(output.codon,'a') as out:
            out.write('>' + str(line[0].description) +'\n')
            out.write(''.join(u))
            out.write('\n')

######################################
### This rule reads in the codon fasta file and also
### the germline amino acid V region fasta files 
### to make a profile alignment with the germline seq
######################################
rule profile_alignment:
  input:
    co=rules.codon_maker.output.codon,
    germ="data/input/Germline_nuc_V{v_gene}.fasta"
  ####
  output:
    profile="data/out/{patient_id}/V{v_gene}_profile.fasta"
  ####
  shell:
    "mafft --add {input.germ} --reorder {input.co} > {output.profile}"  

######################################
### This rule reads in the codon fasta file 
### and trims the gaps in a sequence if there is 
### more than 60% identity at the site
######################################
rule gap_trimmer:
  input:
    cod=rules.codon_maker.output.codon,
  ####
  output:
    ungapped="data/out/{patient_id}/V{v_gene}_ungapped.fasta",
  ####
  run:
    #import pdb;pdb.set_trace()
    import numpy as np
    codon_file = list(SeqIO.parse(input.cod, 'fasta'))
    id_codon = []
    seq_codon = []
    for i in codon_file:
      id_codon.append(str(i.id))
      seq_codon.append(list(i.seq))

    temp = np.array(seq_codon)
    gap_fraction = np.sum(temp=='-', axis=0)/temp.shape[0]
    good_sites = gap_fraction < .4
    temp[:, good_sites]
    list(temp[:, good_sites])
    j = [''.join(row) for row in list(temp[:, good_sites])]

    lines = zip(id_codon,j)
    with open(output.ungapped, 'w' ) as file:
      for line in lines:
        file.write("{}{}\n{}\n".format( '>', line[0], line[1]))

######################################
### This rule builds trees from all the
### ungapped, aligned, codon fasta files
######################################
rule trees:
  input:
    ung=rules.gap_trimmer.output.ungapped
  #### 
  output:
    tree="data/out/{patient_id}/V{v_gene}.new"
  #### 
  # run:
  #   import pdb;pdb.set_trace()
  shell:
    "FastTree -nt {input.ung} > {output.tree}"

######################################
### This rule builds a json for the viz
### which consists of the fasta (which is the profile alignment),
### the newick file (built from the ungapped, aligned, codon fasta files),
### the CDR3 coordinates for each sequence, and the FR3 coordinates
######################################
rule v_gene_json:
  input:
    ali_aa=rules.alignments.output.dog,
    pro=rules.profile_alignment.output.profile,
    tr=rules.trees.output.tree,
    grm=rules.profile_alignment.input.germ
  #### 
  output:
    json="data/out/{patient_id}/V{v_gene}.json"
  #### 
  run:
    from Bio import SeqIO
    #import pdb;pdb.set_trace()
    import numpy as np
    profile = list(SeqIO.parse(input.pro, 'fasta'))
    seq = list(SeqIO.parse(input.grm, 'fasta'))

    ## is this supposed to be a list? 
    profiled_seq = None
    #profiled_seq = []
    for seq_record in profile:
      if 'Germline' in seq_record.description:
        ## and should we be appending to it here??
        profiled_seq = seq_record
        #profiled_seq.append(seq_record)ß
        print('yes')

    profiled_np = np.array(list(str(profiled_seq.seq)), dtype='<U1')
    is_gap = profiled_np == '-'
    profile_indices = np.arange(len(is_gap))
    index_map = profile_indices[~is_gap]

    ## these are nucleotide distances ##
    CDR3_dict = {'V1':(557,564), 'V1-2':(557,564), 'V1-3':(451,458), 'V1-8':(489,496), 'V1-18':(476,483), 'V1-24':(498,505), 'V1-38-4':(679,686), 'V1-45':(432,439), 'V1-46':(583,590), 'V1-58':(581,588), 'V1-68':(6744,6751), 'V1-69':(664,671), 'V1-69-2':(681,688), 'V1-69D':(682,689), 'V2':(505,514), 'V2-5':(505,514), 'V2-10':(502,511), 'V2-26':(455,464), 'V2-70':(435,444), 'V3':(632,639), 'V3-7':(632,639), 'V3-9':(568,577), 'V3-11':(490,497), 'V3-13':(447,454), 'V3-15':(456,463), 'V3-16':(476,483), 'V3-19':(584,591), 'V3-20':(458,465), 'V3-21':(11917,11924), 'V3-22':(539,546), 'V3-23':(458,465), 'V3-25':(524,531), 'V3-29':(11917,11924), 'V3-30':(2228,2235), 'V3-30-3':(2134,2141), 'V3-32':(23127,23136), 'V3-33':(11917,11924), 'V3-35':(586,593), 'V3-38':(451,460), 'V3-43':(504,617), 'V3-43D':(699,708), 'V3-47':(286,291), 'V3-48':(622,629), 'V3-49':(678,685), 'V3-52':(655,662), 'V3-53':(481,488), 'V3-54':(585,592), 'V3-62':(6744,6751), 'V3-63':(458,467), 'V3-64':(529,536), 'V3-66':(445,452), 'V3-69':(454,461), 'V3-71':(6744,6751), 'V3-72':(541,548), 'V3-73':(978,985), 'V3-74':(471,478), 'V4':(580,587), 'V4-4':(580,587), 'V4-28':(578,585), 'V4-30-2':(292,299), 'V4-30-4':(421,438), 'V4-31':(318,325), 'V4-34':(11917,11924), 'V4-38-2':(289,294), 'V4-39':(11917,11924), 'V4-55':(658,665), 'V4-59':(2332,2339), 'V4-61':(581,588), 'V5':(301,306), 'V5-10-1':(301,306), 'V5-51':(596,603), 'V6':(777,784), 'V6-1':(777,784)}
    FR3_dict = {'V1':(443,556),'V1-2':(443,556), 'V1-3':(337,450), 'V1-8':(375,488), 'V1-18':(362,475), 'V1-24':(384,497), 'V1-38-4':(565,678), 'V1-45':(318,431), 'V1-46':(469,582), 'V1-58':(467,580), 'V1-68':(6630,6743), 'V1-69':(550,663), 'V1-69-2':(567,680), 'V1-69D':(568,681), 'V2':(391,504), 'V2-5':(391,504), 'V2-10':(388,501), 'V2-26':(341,454), 'V2-70':(321,434), 'V3':(518,631), 'V3-7':(518,631), 'V3-9':(454,567), 'V3-11':(376,489), 'V3-13':(333,446), 'V3-15':(342,455), 'V3-16':(362,475), 'V3-19':(470,583), 'V3-20':(344,457), 'V3-21':(11803,11916), 'V3-22':(425,538), 'V3-23':(344,457), 'V3-25':(410,523), 'V3-29':(11803,11916), 'V3-30':(2114,2227), 'V3-30-3':(2142,2255), 'V3-32':(23013,23126), 'V3-33':(11803,11916), 'V3-35':(472,585), 'V3-38':(337,450), 'V3-43':(618,627), 'V3-43D':(585,698), 'V3-47':(172,285), 'V3-48':(508,621), 'V3-49':(564,677), 'V3-52':(541,654), 'V3-53':(367,480), 'V3-54':(471,584), 'V3-62':(6630,6743), 'V3-63':(344,457), 'V3-64':(415,528), 'V3-66':(331,444), 'V3-69':(340,453), 'V3-71':(6630,6743), 'V3-72':(427,540), 'V3-73':(864,977), 'V3-74':(357,470), 'V4':(466,579), 'V4-4':(466,579), 'V4-28':(464,577), 'V4-30-2':(178,291), 'V4-30-4':(317,430), 'V4-31':(204,317), 'V4-34':(11803,11916), 'V4-38-2':(175,288), 'V4-39':(11803,11916), 'V4-55':(544,657), 'V4-59':(2218,2331), 'V4-61':(467,580), 'V5':(187,300), 'V5-10-1':(187,300), 'V5-51':(482,595),'V6':(663,776), 'V6-1':(777,784)}

    CDR3_coords = CDR3_dict[str('V'+wildcards.v_gene)]
    CDR3_profile_coords = (int(index_map[CDR3_coords[0]-1]), int(index_map[CDR3_coords[1]-1]))
    FR3_coords = FR3_dict['V'+str(wildcards.v_gene)]
    FR3_profile_coords = (int(index_map[FR3_coords[0]-1]), int(index_map[FR3_coords[1]-1]))

    with open(input.pro) as file:
        fasta = file.read()
    with open(input.tr) as file:
        newick = file.read()

    output_dict = {
        'fasta': fasta,
        'newick': newick,
        'CDR3': CDR3_profile_coords,
        'FR3': FR3_profile_coords
    }

    with open(output.json, 'w') as file:
        json.dump(output_dict, file)

