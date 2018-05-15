# Extract fasta from JSON
python extract_fasta_from_json.py --file data_non_in/7_clone.json --size 30
python extract_fasta_from_json.py --file data_non_in/9_clone.json --size 30
python extract_fasta_from_json.py --file data_non_in/11_clone.json --size 30

cat data_non_in/*clone*_unaligned.fasta > data_non_out/full_size-30_unaligned.fasta
mafft data_non_out/full_size-30_unaligned.fasta > data_non_out/full_size-30_non.fasta
FastTree -nt data_non_out/full_size-30_non.fasta > data_non_out/full_size-30.new

cp data_non_out/full_size-30_non.fasta downstream_data_files/full_size-30_non.fasta

python downstream_data_files/Unique_vdj_non.py --file downstream_data_files/full_size-30_non.fasta

mv seq_non.csv downstream_data_files/seq_non.csv