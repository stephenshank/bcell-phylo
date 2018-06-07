   #align all the v genes at all times and then build a tree for them all   

   #?? do we want to cat all the seqs for a particular time?   

cat V*_7.fasta > all_7.fasta
cat V*_9.fasta > all_9.fasta
cat V*_11.fasta > all_11.fasta

cat V1_7.fasta V1_9.fasta V1_11.fasta > V1_all.fasta
cat V2_7.fasta V2_9.fasta V2_11.fasta > V2_all.fasta
cat V3_7.fasta V3_9.fasta V3_11.fasta > V3_all.fasta
cat V4_7.fasta V4_9.fasta V4_11.fasta > V4_all.fasta
cat V5_7.fasta V5_9.fasta V5_11.fasta > V5_all.fasta
cat V6_9.fasta V6_11.fasta > V6_all.fasta

for i in {1..6}
do
  mafft V$i\_all.fasta > V$i\_all_aligned.fasta
done;

for i in {1..6}
do
  FastTree -nt V$i\_all_aligned.fasta > 77612_V$i\_all.new
done;

mafft all_7.fasta > 77612_all_7_aligned.fasta
mafft all_9.fasta > 77612_all_9_aligned.fasta
mafft all_11.fasta > 77612_all_11_aligned.fasta

mafft V1_7.fasta > V1_7_aligned.fasta
mafft V1_9.fasta > V1_9_aligned.fasta
mafft V1_11.fasta > V1_11_aligned.fasta
mafft V2_7.fasta > V2_7_aligned.fasta
mafft V2_9.fasta > V2_9_aligned.fasta
mafft V2_11.fasta > V2_11_aligned.fasta
mafft V3_7.fasta > V3_7_aligned.fasta
mafft V3_9.fasta > V3_9_aligned.fasta
mafft V3_11.fasta > V3_11_aligned.fasta
mafft V4_7.fasta > V4_7_aligned.fasta
mafft V4_9.fasta > V4_9_aligned.fasta
mafft V4_11.fasta > V4_11_aligned.fasta
mafft V5_7.fasta > V5_7_aligned.fasta
mafft V5_9.fasta > V5_9_aligned.fasta
mafft V5_11.fasta > V5_11_aligned.fasta
mafft V6_7.fasta > V6_9_aligned.fasta
mafft V6_9.fasta > V6_9_aligned.fasta
mafft V6_11.fasta > V6_11_aligned.fasta

FastTree -nt V1_7_aligned.fasta > 77612_V1_7.new
FastTree -nt V1_9_aligned.fasta > 77612_V1_9.new
FastTree -nt V1_11_aligned.fasta > 77612_V1_11.new
FastTree -nt V2_7_aligned.fasta > 77612_V2_7.new
FastTree -nt V2_9_aligned.fasta > 77612_V2_9.new
FastTree -nt V2_11_aligned.fasta > 77612_V2_11.new
FastTree -nt V3_7_aligned.fasta > 77612_V3_7.new
FastTree -nt V3_9_aligned.fasta > 77612_V3_9.new
FastTree -nt V3_11_aligned.fasta > 77612_V3_11.new
FastTree -nt V4_7_aligned.fasta > 77612_V4_7.new
FastTree -nt V4_9_aligned.fasta > 77612_V4_9.new
FastTree -nt V4_11_aligned.fasta > 77612_V4_11.new
FastTree -nt V5_7_aligned.fasta > 77612_V5_7.new
FastTree -nt V5_9_aligned.fasta > 77612_V5_9.new
FastTree -nt V5_11_aligned.fasta > 77612_V5_11.new
FastTree -nt V6_7_aligned.fasta > 77612_V6_7.new
FastTree -nt V6_9_aligned.fasta > 77612_V6_9.new
FastTree -nt V6_11_aligned.fasta > 77612_V6_11.new

FastTree -nt 77612_all_7_aligned.fasta > 77612_all_7.new
FastTree -nt 77612_all_9_aligned.fasta > 77612_all_9.new
FastTree -nt 77612_all_11_aligned.fasta > 77612_all_11.new
