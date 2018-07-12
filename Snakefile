PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
GENES = ["V1", "V2", "V3", "V4", "V5", "V6", "V7"] 

rule all:
  input:
    expand("data/out/{patient_id}/{v_gene}.json", patient_id=PATIENT_IDS, v_gene=GENES)

rule clone_unaligned_fasta:
  input:
    "data/input/{patient_id}_{clone}_clone.json",
    "python/json_to_unaligned_fasta.py"
  output:
    "data/out/{patient_id}/clone_{clone}_unaligned.fasta"
  shell:
    "python python/json_to_unaligned_fasta.py -p {wildcards.patient_id} -c {wildcards.clone} --size 30"

rule gene_unaligned_fasta:
  input:
    "data/out/{patient_id}/clone_1_unaligned.fasta",
    "data/out/{patient_id}/clone_2_unaligned.fasta",
    "data/out/{patient_id}/clone_3_unaligned.fasta",
    "data/out/{patient_id}/clone_4_unaligned.fasta",
    "data/out/{patient_id}/clone_5_unaligned.fasta",
    "data/out/{patient_id}/clone_6_unaligned.fasta",
    "python/full_fasta_to_v-gene_separate_fasta.py"
  output:
    "data/out/{patient_id}/V1_unaligned.fasta",
    "data/out/{patient_id}/V2_unaligned.fasta",
    "data/out/{patient_id}/V3_unaligned.fasta",
    "data/out/{patient_id}/V4_unaligned.fasta",
    "data/out/{patient_id}/V5_unaligned.fasta",
    "data/out/{patient_id}/V6_unaligned.fasta",
    "data/out/{patient_id}/V7_unaligned.fasta"
  shell:
    "python python/full_fasta_to_v-gene_separate_fasta.py -p {wildcards.patient_id}"

rule alignments:
  input:
    "data/out/{patient_id}/V{v_gene}_unaligned.fasta"
  output:
    "data/out/{patient_id}/V{v_gene}.fasta"
  shell:
    "mafft {input} > {output}"

rule trees:
  input:
    "data/out/{patient_id}/V{v_gene}.fasta"
  output:
    "data/out/{patient_id}/V{v_gene}.new"
  shell:
    "FastTree -nt {input} > {output}"

rule v_gene_json:
  input:
    "data/out/{patient_id}/V{v_gene}.fasta",
    "data/out/{patient_id}/V{v_gene}.new",
    "python/json_for_dashboard.py"
  output:
    "data/out/{patient_id}/V{v_gene}.json"
  shell:
    "python python/json_for_dashboard.py -p {wildcards.patient_id} -g {wildcards.v_gene}"
