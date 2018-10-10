from python import *


PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
CLONES = ["1", "2", "3", "4", "5", "6"]
GENES = [
  "1-18", "1-2", "1-24", "1-3", "1-46", "1-69", "1-69-2", "1-69D", "1-8", "1-45", "1-58",
  "2-5", "2-70", "2-26",
  "3-11", "3-13", "3-15", "3-20", "3-21", "3-23", "3-30-3", "3-33", "3-43", "3-48", "3-49", "3-53", "3-64", "3-66", "3-7", "3-72", "3-9", "3-43D", "3-73", "3-74",
  "4-30-2", "4-30-4", "4-31", "4-34", "4-38-2", "4-39", "4-4", "4-59", "4-61",
  "5-51", "5-10-1",
  "6-1"
]

rule all:
  input:
    expand(
      "data/{patient_id}/V{v_gene}.json",
      patient_id=PATIENT_IDS,
      v_gene=GENES,
      clone=CLONES
    )

rule unpacked:
  input:
    "data/bcell-phylo_Ver3.tar.gz"
  output:
    expand("data/input/{patient_id}_{clone}_clone.json", patient_id=PATIENT_IDS, clone=CLONES),
    expand("data/input/Germline_nuc_V{gene}.fasta", gene=GENES)
  shell:
    "tar xvzf {input} -C data/input"

rule clone_json_to_unaligned_fasta:
  input:
    json="data/input/{patient_id}_{clone}_clone.json"
  output:
    fasta="data/{patient_id}/clone_{clone}_unaligned.fasta"
  run:
    clone_json_to_unaligned_fasta(input.json, output.fasta, wildcards.clone)

rule separate_into_regions:
  input:
    fasta=expand("data/{{patient_id}}/clone_{clone}_unaligned.fasta", clone=CLONES)
  output:
    fasta="data/{patient_id}/V{v_gene}_unaligned.fasta"
  run:
    separate_into_regions(input.fasta, output.fasta, wildcards.v_gene)

rule protein_and_corrected_dna:
  input:
    fasta="data/{patient_id}/V{v_gene}_unaligned.fasta",
  output:
    aa="data/{patient_id}/V{v_gene}_unaligned_corrected_AA.fasta",
    nuc="data/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta"
  run:
    protein_and_corrected_dna(input.fasta, output.aa, output.nuc)

rule alignments:
  input:
    rules.protein_and_corrected_dna.output.aa
  output:
    "data/{patient_id}/V{v_gene}_aligned_AA.fasta",
  shell:
    "mafft --amino {input} > {output}"

rule codon_maker:
  input:
    aligned_aa=rules.alignments.output[0],
    unaligned_nucleotide=rules.protein_and_corrected_dna.output.nuc
  output:
    codon="data/{patient_id}/V{v_gene}_codon.fasta",
  run:
    protein_alignment_to_codon_alignment(input.aligned_aa, input.unaligned_nucleotide, output.codon)

rule profile_alignment:
  input:
    codon=rules.codon_maker.output.codon,
    germline="data/input/Germline_nuc_V{v_gene}.fasta"
  output:
    profile="data/{patient_id}/V{v_gene}_profile.fasta"
  shell:
    "mafft --add {input.germline} --reorder {input.codon} > {output.profile}"  

rule gap_trimmer:
  input:
    codon=rules.codon_maker.output.codon,
  output:
    trimmed="data/{patient_id}/V{v_gene}_ungapped.fasta",
  run:
    gap_trimmer(input.codon, output.trimmed)

rule trees:
  input:
    fasta=rules.gap_trimmer.output.trimmed
  output:
    tree="data/{patient_id}/V{v_gene}.new"
  shell:
    "FastTree -nt {input.fasta} > {output.tree}"

rule v_gene_json:
  input:
    profile=rules.profile_alignment.output.profile,
    tree=rules.trees.output.tree,
    germline=rules.profile_alignment.input.germline
  output:
    json="data/{patient_id}/V{v_gene}.json"
  run:
    json_for_dashboard(input.profile, input.tree, input.germline, output.json, wildcards)

