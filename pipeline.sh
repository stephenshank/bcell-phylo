# Extract fasta from JSON
python extract_fasta_from_json.py --file data/6_clone.json --size 30
python extract_fasta_from_json.py --file data/7_clone.json --size 30
python extract_fasta_from_json.py --file data/9_clone.json --size 30
python extract_fasta_from_json.py --file data/11_clone.json --size 30

cat data/*clone*_unaligned.fasta > data/full_size-30_unaligned.fasta
mafft data/full_size-30_unaligned.fasta > data/full_size-30.fasta
FastTree -nt data/full_size-30.fasta > data/full_size-30.new
