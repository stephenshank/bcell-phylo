from python import *


PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
CLONES = ["1", "2", "3", "4", "5", "6"]
GENES = [
  "1", "1-18", "1-2", "1-24", "1-3", "1-46", "1-69", "1-69-2", "1-69D", "1-8", "1-45", "1-58",
  "2", "2-5", "2-70", "2-26",
  "3", "3-11", "3-13", "3-15", "3-20", "3-21", "3-23", "3-30-3", "3-33", "3-43", "3-48", "3-49", "3-53", "3-64", "3-66", "3-7", "3-72", "3-9", "3-43D", "3-73", "3-74",
  "4", "4-30-2", "4-30-4", "4-31", "4-34", "4-38-2", "4-39", "4-4", "4-59", "4-61",
  "5", "5-51", "5-10-1",
  "6", "6-1"
]
IMGT_IDS = [
  "M99637", "M99641", "M99642", "M99652",  "M99653",  "M99654"
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
    "data/bcell-phylo_Ver4.tar.gz"
  output:
    expand("data/input/{patient_id}_{clone}_clone.json", patient_id=PATIENT_IDS, clone=CLONES),
    "data/input/imgt.dat",
    "data/input/imgt.fasta"
  shell:
    "tar xvzf {input} -C data/input"

rule clone_json_to_unaligned_fasta:
  input:
    json="data/input/{patient_id}_{clone}_clone.json"
  output:
    fasta="data/{patient_id}/clone_{clone}_unaligned.fasta"
  run:
    clone_json_to_unaligned_fasta(input.json, output.fasta, wildcards.clone)

rule unique_vs:
  input:
    expand("data/input/{patient_id}_{clone}_clone.json", patient_id=PATIENT_IDS, clone=CLONES),
  output:
    "data/unique_vs.json"
  run:
    get_unique_vs(PATIENT_IDS, CLONES)

rule imgt_records:
  input:
    dat="data/input/imgt.dat"
  output:
    "data/imgt/V{v_gene}/raw.txt"
  run:
    extract_imgt_records(input.dat, IMGT_IDS)

rule imgt_information:
  input:
    raw="data/imgt/V{v_gene}/raw.txt"
  output:
    nucleotide_fasta="data/imgt/V{v_gene}/nucleotide.fasta",
    protein_fasta="data/imgt/V{v_gene}/protein.fasta",
    json="data/imgt/V{v_gene}/data.json"
  run:
    parse_imgt_record(input.raw, output.nucleotide_fasta, output.protein_fasta, output.json, wildcards.v_gene)

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
    codon=rules.alignments.output[0],
    germline=rules.imgt_information.output.protein_fasta
  output:
    fasta="data/{patient_id}/V{v_gene}_profile.fasta"
  shell:
    "mafft --add {input.germline} --reorder {input.codon} > {output.fasta}"  

rule gap_trimmer:
  input:
    codon=rules.codon_maker.output.codon
  output:
    trimmed="data/{patient_id}/V{v_gene}_ungapped.fasta",
  run:
    gap_trimmer(input.codon, output.trimmed)

rule indicial_mapper:
  input:
    fasta=rules.profile_alignment.output.fasta,
    json=rules.imgt_information.output.json
  output:
    json="data/{patient_id}/V{v_gene}_indices.json"
  run:
    indicial_mapper(input.fasta, input.json, output.json, wildcards.v_gene)

rule trees:
  input:
    fasta=rules.gap_trimmer.output.trimmed
  output:
    tree="data/{patient_id}/V{v_gene}.new"
  shell:
    "FastTree -nt {input.fasta} > {output.tree}"

rule v_gene_json:
  input:
    fasta=rules.profile_alignment.output.fasta,
    json=rules.indicial_mapper.output.json,
    tree=rules.trees.output.tree,
  output:
    json="data/{patient_id}/V{v_gene}.json"
  run:
    json_for_dashboard(input.fasta, input.json, input.tree, output.json, wildcards)

