# Extract fasta from JSON
for i in {7..12}
do
  python python/extract_fasta_from_json.py \
    --file data/input/$i\_clone.json \
    --size 30 \
    --output data/replicates;
done;

cat data/replicates/*clone*_unaligned.fasta > data/replicates/full_size-30_unaligned.fasta
mafft data/replicates/full_size-30_unaligned.fasta > data/replicates/full_size-30.fasta
FastTree -nt data/replicates/full_size-30.fasta > data/replicates/full_size-30.new

#python python/extract_csv_from_fasta.py --file data/replicates/full_size-30.fasta
