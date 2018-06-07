## this is just for one json ##

# use as: bash pipeline/v1.sh
#python python/try1.py --file data/9_clone_size-30_unaligned.fasta --output data/v1

#mafft data/v1/V1_un.fasta > data/v1/V1_aligned.fasta

#######################################
## THIS IS WHERE THE PIPELINE STARTS ##
######################################

# use as: bash pipeline/v_get.sh

	## this loop will grab the jsons and spit out the clone files (with all the sequences) ##
	## need to change this loop to match the patient json files ## 

## this is an example of a start stop step loop in bash ##
#for i in $(seq 4 2 8)
#do
#   echo "Welcome $i times"
#done

#for i in $(seq 4 2 8)
for i in {7..12}
do
  python python/vfull_json_to_fasta.py \
	--file data/input_json/$i\_clone.json \
 	--size 30 \
	--output data/output_unaligned_clone_files;
done;

	## this will cat all the unaligned clone files together into one 'full' file ## 
cat data/output_unaligned_clone_files/*_clone_size-30_unaligned.fasta > data/output_unaligned_clone_files/full_clone_size-30_unaligned.fasta

	## this loop will take the unaligned clone files and take all the different v genes and make unaligned files for each gene ##
	## need to change this loop to match the patient json files ##
#for i in $(seq 4 2 8)
for i in {7..12}
do
  python python/full_fasta_to_v-gene_separate_fasta.py \
	--file data/output_unaligned_clone_files/$i\_clone_size-30_unaligned.fasta \
	--output data/output_unaligned_v-genes;
done;

	## this is aligning the 'full' clone file ##
mafft data/output_unaligned_clone_files/full_clone_size-30_unaligned.fasta > data/output_aligned_clone_files/full_aligned.fasta


	## this is aligning the first V gene that (does not get assigned?)## 
mafft data/output_unaligned_v-genes/V_un.fasta > data/output_aligned_v-genes/V_aligned.fasta

	## this loop grabs all the unaligned gene files and then aligns them ##
for i in {1..7}
do
  mafft data/output_unaligned_v-genes/V$i\_un.fasta > data/output_aligned_v-genes/V$i\_aligned.fasta

done;

	## this will take the 'full' aligned fasta file, representing sequences from all the time points and make a tsv ##
python python/full_aligned_to_tsv.py \
	--file data/output_aligned_clone_files/full_aligned.fasta \
	--output data/output_aligned_clone_files

	## this builds a tree for the 'full' aligned file ##
FastTree -nt data/output_aligned_clone_files/full_aligned.fasta > data/output_newick/77612_full.new
	
	## this builds a tree for the unassigned v gene sequences ##
FastTree -nt data/output_aligned_v-genes/V_aligned.fasta > data/output_newick/77612_V.new

	## this loop builds a newick tree for all the other gene sequences that were collected ## 
for i in {1..7}
do
  FastTree -nt data/output_aligned_v-genes/V$i\_aligned.fasta > data/output_newick/77612_V$i\.new;
done;


Rscript R/R_viz_patient.r data/output_aligned_clone_files/77612_master.tsv data/pretty_pictures

