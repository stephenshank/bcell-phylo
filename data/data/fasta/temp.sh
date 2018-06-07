mafft V6_7.fasta > V6_7_aligned.fasta
mafft V6_9.fasta > V6_9_aligned.fasta
mafft V6_11.fasta > V6_11_aligned.fasta

FastTree -nt V6_7_aligned.fasta > 77612_V6_7.new
FastTree -nt V6_9_aligned.fasta > 77612_V6_9.new
FastTree -nt V6_11_aligned.fasta > 77612_V6_11.new