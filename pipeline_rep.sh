# Extract fasta from JSON
python extract_fasta_from_json.py --file data_rep_in/7_clone.json --size 30
python extract_fasta_from_json.py --file data_rep_in/9_clone.json --size 30
python extract_fasta_from_json.py --file data_rep_in/11_clone.json --size 30
python extract_fasta_from_json.py --file data_rep_in/8_clone.json --size 30
python extract_fasta_from_json.py --file data_rep_in/10_clone.json --size 30
python extract_fasta_from_json.py --file data_rep_in/12_clone.json --size 30

# the last thing I want to do is check
cat data_rep_in/*clone*_unaligned.fasta > data_rep_out/full_size-30_unaligned.fasta
mafft data_rep_out/full_size-30_unaligned.fasta > data_rep_out/full_size-30_rep.fasta
FastTree -nt data_rep_out/full_size-30_rep.fasta > data_rep_out/full_size-30.new

cp data_rep_out/full_size-30_rep.fasta downstream_data_files/full_size-30_rep.fasta

python downstream_data_files/Unique_vdj_rep.py --file downstream_data_files/full_size-30_rep.fasta

mv seq_rep.csv downstream_data_files/seq_rep.csv