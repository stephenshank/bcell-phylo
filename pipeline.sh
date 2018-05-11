# Extract fasta from JSON
python extract_fasta_from_json.py --file data/1_clone.json --size 30
python extract_fasta_from_json.py --file data/9_clone.json --size 30
python extract_fasta_from_json.py --file data/18_clone.json --size 30
python extract_fasta_from_json.py --file data/27_clone.json --size 30

cat data/*clone*_unaligned.fasta > data/full_size-30_unaligned.fasta
mafft data/full_size-30_unaligned.fasta > data/full_size-30.fasta
FastTree -nt data/full_size-30.fasta > data/full_size-30.new
