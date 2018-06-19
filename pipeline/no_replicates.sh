# Extract fasta from JSON
for size in 30 200; do
  python python/extract_fasta_from_json.py --file data/input/7_clone.json --size $size --output data/no_replicates
  python python/extract_fasta_from_json.py --file data/input/9_clone.json --size $size --output data/no_replicates
  python python/extract_fasta_from_json.py --file data/input/11_clone.json --size $size --output data/no_replicates

  cat data/no_replicates/*clone*_size-$size\_unaligned.fasta > data/no_replicates/full_size-$size\_unaligned.fasta
  remove_duplicates data/no_replicates/full_size-$size\_unaligned.fasta
  mafft data/no_replicates/full_size-$size\_unaligned-no_duplicates.fasta > data/no_replicates/full_size-$size\.fasta
  FastTree -nt data/no_replicates/full_size-$size\.fasta > data/no_replicates/full_size-$size\.new

  python python/extract_csv_from_fasta.py --file data/no_replicates/full_size-$size\.fasta
done
