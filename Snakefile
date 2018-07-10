PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
GENES = ["V1", "V2", "V3", "V4", "V5", "V6", "V7"] 

rule all:
  input:
    expand("data/out/{patient_id}/{gene}.new", patient_id=PATIENT_IDS, gene=GENES)

rule unaligned_fasta_clone:
  input:
    "data/input/{patient_id}_{clone}_clone.json"
  output:
    "data/out/{patient_id}/{patient_id}_{clone}_clone_size-30_unaligned.fasta"
  shell:
    "python python/json_to_unaligned_fasta.py --file {input} --size 30 --out data/out/{wildcards.patient_id}"

