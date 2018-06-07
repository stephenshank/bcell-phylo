python ../python/full_aligned_to_tsv.py \
	--file output_aligned_clone_files/full_aligned.fasta \
	--output output_aligned_clone_files/
	
Rscript ../R/R_viz_patient.r 77612_master.tsv pretty_pictures/
