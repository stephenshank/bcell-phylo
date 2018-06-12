	## this builds a tree for the 'full' aligned file ##
FastTree -nt data/sandbox/output/aligned_jsons/full_aligned.fasta > data/sandbox/viz/trees/77612_full.new
	
	## this builds a tree for the unassigned v gene sequences ##
FastTree -nt data/sandbox/output/aligned_v-genes/V_aligned.fasta > data/sandbox/viz/trees/77612_V.new

	## this loop builds a newick tree for all the other gene sequences that were collected ## 
for i in {1..7}
do
  FastTree -nt data/sandbox/output/aligned_v-genes/V$i\_aligned.fasta > data/sandbox/viz/trees/77612_V$i\.new;
done;

	## this loop will align all the time-point aligned fasta files ##
for i in {7..12}
do
  FastTree -nt data/sandbox/output/aligned_jsons/$i\_aligned.fasta > data/sandbox/viz/trees/77612_$i\.new;
done;