# Extract fasta from JSON
python python/extract_fasta_from_json.py --file data/input/7_clone.json --size 30 --output data/no_replicates
python python/extract_fasta_from_json.py --file data/input/9_clone.json --size 30 --output data/no_replicates
python python/extract_fasta_from_json.py --file data/input/11_clone.json --size 30 --output data/no_replicates

cat data/no_replicates/*clone*_unaligned.fasta > data/no_replicates/full_size-30_unaligned.fasta
mafft data/no_replicates/full_size-30_unaligned.fasta > data/no_replicates/full_size-30.fasta
FastTree -nt data/no_replicates/full_size-30.fasta > data/no_replicates/full_size-30.new

python python/extract_csv_from_fasta.py --file data/no_replicates/full_size-30.fasta
