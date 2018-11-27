from python import *


PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
CLONES = ["1", "2", "3", "4", "5", "6"]
UNIQUE_VS, PATIENT_V_PAIRS = get_unique_vs(PATIENT_IDS, CLONES)

"""
rule all:
  input:
    expand(
      "data/{patient_id}/V{v_gene}.json",
      patient_id=PATIENT_IDS,
      v_gene=GENES,
    )
"""

rule imgt_information:
  input:
    data="data/input/imgt.dat"
  output:
    fasta="data/imgt/sequences.fasta"
  run:
    extract_imgt_records(input.data)
    get_imgt_information(input.data, output.fasta)

rule blast_database:
  input:
    rules.imgt_information.output.fasta
  output:
    "data/blast"
  shell:
    "makeblastdb -in {input} -dbtype nucl -parse_seqids -out data/blast/db"

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
    fasta="data/{patient_id}/V{v_gene}/unaligned.fasta"
  run:
    separate_into_regions(input.fasta, output.fasta, wildcards.v_gene)

rule collapse_identical_sequences:
  input:
    fasta=rules.separate_into_regions.output.fasta,
  output:
    fasta="data/{patient_id}/V{v_gene}/unaligned_collapsed.fasta"
  run:
    collapse_identical_sequences(input.fasta, output.fasta)

rule protein_and_corrected_dna:
  input:
    fasta=rules.collapse_identical_sequences.output.fasta
  output:
    aa="data/{patient_id}/V{v_gene}/unaligned_corrected_AA.fasta",
    nuc="data/{patient_id}/V{v_gene}/unaligned_corrected_nuc.fasta"
  run:
    protein_and_corrected_dna(input.fasta, output.aa, output.nuc)

rule alignments:
  input:
    rules.protein_and_corrected_dna.output.aa
  output:
    "data/{patient_id}/V{v_gene}/aligned_AA.fasta",
  shell:
    "mafft --amino {input} > {output}"

rule codon_maker:
  input:
    aligned_aa=rules.alignments.output[0],
    unaligned_nucleotide=rules.protein_and_corrected_dna.output.nuc
  output:
    codon="data/{patient_id}/V{v_gene}/codon.fasta",
  run:
    protein_alignment_to_codon_alignment(input.aligned_aa, input.unaligned_nucleotide, output.codon)

rule consensus:
  input:
    fasta=rules.codon_maker.output.codon
  output:
    fasta="data/{patient_id}/V{v_gene}/consensus.fasta"
  run:
    get_consensus_sequence(input.fasta, output.fasta)

rule blast_results:
  input:
    rules.consensus.output.fasta
  output:
    "data/{patient_id}/V{v_gene}/blast/result.json"
  shell:
    "blastn -db data/blast/db -outfmt 13 -query {input} -out {output}"

"""
rule profile_alignment:
  input:
    codon=rules.alignments.output[0],
    germline=rules.imgt_information.output.protein_fasta
  output:
    fasta="data/{patient_id}/V{v_gene}_profile.fasta"
  shell:
    "mafft --add {input.germline} --reorder {input.codon} > {output.fasta}"  
"""
rule gap_trimmer:
  input:
    codon=rules.codon_maker.output.codon
  output:
    trimmed="data/{patient_id}/V{v_gene}/codon_ungapped.fasta",
  run:
    gap_trimmer(input.codon, output.trimmed)
"""
rule indicial_mapper:
  input:
    fasta=rules.profile_alignment.output.fasta,
    json=rules.imgt_information.output.json
  output:
    json="data/{patient_id}/V{v_gene}_indices.json"
  run:
    indicial_mapper(input.fasta, input.json, output.json, wildcards.v_gene)
"""
rule trees:
  input:
    fasta=rules.gap_trimmer.output.trimmed
  output:
    tree="data/{patient_id}/V{v_gene}/tree.new"
  shell:
    "FastTree -nt {input.fasta} > {output.tree}"

"""
rule v_gene_json:
  input:
    fasta=rules.profile_alignment.output.fasta,
    json=rules.indicial_mapper.output.json,
    tree=rules.trees.output.tree,
  output:
    json="data/{patient_id}/V{v_gene}/dashboard.json"
  run:
    json_for_dashboard(input.fasta, input.json, input.tree, output.json, wildcards)
"""
