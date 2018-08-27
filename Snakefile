include: "rules/common.smk"

#PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
#GENES = ["1", "2", "3", "4", "5", "6"] 
#

def debug(wildcards):
  import pdb;pdb.set_trace()

# TODO: Place in config
PATIENT_IDS = ["28729"]
CLONES = ["1"]
GENES = ["1", "2", "3", "4", "5", "6"]

rule all:
  input:
    #expand("data/out/{patient_id}/{v_gene}_unaligned.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
    #expand("data/out/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)
    expand("data/out/{patient_id}/V{v_gene}_aligned_AA.fasta", patient_id=PATIENT_IDS, v_gene=GENES, clone=CLONES)

rule clone_unaligned_fasta:
  input:
    "data/input/{patient_id}_{v_gene}_clone.json",
    "python/json_to_unaligned_fasta.py"
  output:
    "data/out/{patient_id}/clone_{v_gene}_unaligned.fasta"
  shell:
    "python python/json_to_unaligned_fasta.py -p {wildcards.patient_id} -c {wildcards.v_gene} --size 30"

rule gene_unaligned_fasta:
# We actually want to include all clone files.
  input:
    clones=get_clones
  output:
    v_gene_output="data/out/{patient_id}/V{v_gene}_unaligned.fasta"
  log:
    "logs/{patient_id}_{v_gene}_clone_to_gene.log"
  run:
    v_gene_regex = re.compile('._V(' + wildcards.v_gene + ').')
    bad_searches = []
    bad_search_text = ""
    v_sequences = []
    for f in input.clones:
      current_clone = SeqIO.parse(f, 'fasta')
      for sequence in current_clone:
          regex_search = v_gene_regex.search(sequence.id)
          if not regex_search is None:
              v_gene = int(regex_search.group(1))
              v_sequences.append(sequence)
          else:
              bad_searches.append(sequence.id)
          SeqIO.write(v_sequences, output.v_gene_output, 'fasta')
          print('%d searches were bad: %s' % ( len(bad_searches), bad_search_text))

### again, check the wild cards here ####
rule unaligned_to_amino_acids:
  input:
    "data/out/{patient_id}/V{v_gene}_unaligned.fasta",
    "python/full_gene_trans_to_AA.py"
  output:
    unaligned_aa="data/out/{patient_id}/V{v_gene}_unaligned_corrected_AA.fasta",
    unaligned_nucleotide="data/out/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta"
  shell:
    "python python/full_gene_trans_to_AA.py -p {wildcards.patient_id} -g {wildcards.v_gene}"

############## testing this out with python inside
######################## probably wont work now 
rule alignments:
  input:
    rules.unaligned_to_amino_acids.output.unaligned_aa
  output:
    "data/out/{patient_id}/V{v_gene}_aligned_AA.fasta",
  shell:
    "mafft --amino {input} > {output}"

#rule codon_maker:
#  input:
#    "data/out/{patient_id}/V{v_gene}_aligned_AA.fasta",
#    "data/out/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta",
#    "python/AA_to_codon.py"
#  output:
#    "data/out/{patient_id}/V{v_gene}_codon.fasta",
#  shell:
#    "python python/AA_to_codon.py -p {wildcards.patient_id} -g {wildcards.v_gene}"
#    
#rule profile_alignment:
#  input:
#    "data/out/{patient_id}/V{v_gene}_codon.fasta",
#    "data/input/Germline_nuc_V{v_gene}.fasta"
#  output:
#    "data/out/{patient_id}/V{v_gene}_profile.fasta"
#  shell:
#    "mafft --add data/input/Germline_nuc_V{wildcards.v_gene}.fasta --reorder {input[0]} > {output}"
#    
#rule gap_trimmer:
#  input:
#    "data/out/{patient_id}/V{v_gene}_codon.fasta",
#    "python/gap_trimmer.py"
#  output:
#    "data/out/{patient_id}/V{v_gene}_ungapped.fasta",
#  shell:
#    "python python/gap_trimmer.py -p {wildcards.patient_id} -g {wildcards.v_gene}"
#
#rule trees:
#  input:
#    "data/out/{patient_id}/V{v_gene}_ungapped.fasta"
#  output:
#    "data/out/{patient_id}/V{v_gene}.new"
#  shell:
#    "FastTree -nt {input} > {output}"
#
#rule v_gene_json:
#  input:
#    "data/out/{patient_id}/V{v_gene}_AA.fasta",
#    "data/out/{patient_id}/V{v_gene}_profile.fasta",
#    "data/out/{patient_id}/V{v_gene}.new",
#    "python/json_for_dashboard.py"
#  output:
#    "data/out/{patient_id}/V{v_gene}.json"
#  shell:
#    "python python/json_for_dashboard.py -p {wildcards.patient_id} -g {wildcards.v_gene}"
