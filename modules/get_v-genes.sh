

## this loop will take the unaligned clone files and take all the different v genes and make unaligned files for each gene ##
	## need to change this loop to match the patient json files ##
#for i in $(seq 4 2 8)
for i in {7..12}
do
  python python/full_fasta_to_v-gene_separate_fasta.py \
	--file data/sandbox/output/unaligned_jsons/$i\_clone_size-30_unaligned.fasta \
	--output data/sandbox/output/unaligned_v-genes;
done;


