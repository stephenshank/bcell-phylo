## this is the full list, we will just use the first two to try ##
ids=('28729' '48689' '67029' '77612' '78202' '93954' '99361' '99682' 'GJS')
#ids=('28729' '48689')
#echo ${ids[@]}


######################################
## this first loop will go through each patient ##
for i in ${ids[@]}
do 	
	mkdir data/out/$i
	## this remover will get rid of any files in patients folder ## 
	echo 'removing populated folders for' $i
	rm data/out/$i/$i*.tsv
	rm data/out/$i/$i*.fasta
	rm data/out/$i/$i*.new
 	echo $i 'is a patient'
----####################################
	## this will go through each time point for each patient ## 
	for a in {1..6}
	do
		echo 'converting'$i'_'$a'_clone.json to an unaligned fasta file'
	#######	change the output of this, the size is too rigid ## 
		python python/json_to_unaligned_fasta.py \
				--file data/input/*/$i\_$a\_clone.json \
 				--size 100 \
				--out data/out/$i/
		## this file should be open to append to, but I am curious to see how it 
		## does with multiple input/test-clones/ files..## 
		echo $i $a 'geting v genes from fasta files'	
	###### need to check this python script and the out/$i/put is rigid ## 
		python python/full_fasta_to_v-gene_separate_fasta.py \
				--file data/out/$i/$i\_$a\_*_unaligned.fasta \
				--out data/out/$i/
			
	done;
	echo 'cat all 6 time points'
  ####### cat ing all the fastas made from the jsons at each time point (1-6) ## 
	cat data/out/$i/$i\_*_*_*_unaligned.fasta > data/out/$i/$i\_all_unaligned.fasta
	echo 'making a nice tsv'
	
 #### can use the unaligned fasta file
	python python/full_unaligned_tsv.py \
		--file data/out/$i/$i\_all_unaligned.fasta \
		--out data/out/$i/
	
	echo 'making ggplot of the V gene usage over time points'
	## use the first .split from the file input ## 
	Rscript R/R_viz_patient.r data/out/$i/$i\_master.tsv data/out/pretty_pictures/
	##############################################
	## this will align the gene and build the tree 
	for j in {1..7}
	do
		echo $i $j 'aligning gene'
		mafft data/out/$i/V$j\_un.fasta > data/out/$i/V$j\_aligned.fasta
		echo $i $j 'building tree with this aligned gene'
		FastTree -nt data/out/$i/V$j\_aligned.fasta > data/out/$i/$i\_V$j.new
	done;
done;	
	