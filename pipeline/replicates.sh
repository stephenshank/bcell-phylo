# Extract fasta from JSON
for size in 30 200; do
  for i in {7..12}
  do
    python python/extract_fasta_from_json.py \
      --file data/input/$i\_clone.json \
      --size $size \
      --output data/replicates;
  done;

  cat data/replicates/*clone*_size-$size\_unaligned.fasta > data/replicates/full_size-$size\_unaligned.fasta
  remove_duplicates data/replicates/full_size-$size\_unaligned.fasta
  mafft data/replicates/full_size-$size\_unaligned-no_duplicates.fasta > data/replicates/full_size-$size.fasta
  FastTree -nt data/replicates/full_size-$size.fasta > data/replicates/full_size-$size.new

  python python/extract_csv_from_fasta.py --file data/replicates/full_size-$size.fasta
done
