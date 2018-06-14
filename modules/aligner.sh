
## this is aligning the 'full' clone file ##
## not sure if i even need this alignment... ##
mafft data/sandbox/output/unaligned_jsons/full_clone_size-30_unaligned.fasta > data/sandbox/output/aligned_jsons/full_aligned.fasta

	## this loop is going to align all the time-point fasta files ##
#for i in {7..12}
#do
#  mafft data/sandbox/output/unaligned_jsons/$i\_clone_size-30_unaligned.fasta > data/sandbox/output/aligned_jsons/$i\_aligned.fasta

#done;

	## this is aligning the first V gene that (does not get assigned?)## 
mafft data/sandbox/output/unaligned_v-genes/V_un.fasta > data/sandbox/output/aligned_v-genes/V_aligned.fasta

	## this loop grabs all the unaligned gene files and then aligns them ##
#for i in {1..7}
#do
mafft data/sandbox/output/unaligned_v-genes/V$i\_un.fasta > data/sandbox/output/aligned_v-genes/V$i\_aligned.fasta

#done;
