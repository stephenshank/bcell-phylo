from Bio import SeqIO

def get_clones(wildcards):
  return ["data/out/" + wildcards.patient_id + "/clone_" + clone + "_unaligned.fasta" for clone in CLONES]

def get_rearrangements(wildcards):
  # read unique rearrangements file for each gene
  directory =  "data/out/" + wildcards.patient_id + "/"
  gene_fn =  directory + wildcards.v_gene + "_unique.txt"
  rearrangements_prefixes = open(gene_fn).readlines()
  rearrangement_fns = [directory + prefix.strip() +  "_unaligned_corrected_AA.fasta" for prefix in rearrangements_prefixes]
  return rearrangement_fns


