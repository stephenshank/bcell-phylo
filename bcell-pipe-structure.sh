ids=('28729' '48689' '67029' '77612' '78202' '93954' '99361' '99682' 'GJS')

#echo ${ids[@]}


######################################
## this first loop will go through each patient, so what scripts need this? ##
for i in ${ids[@]}
do 	
	## this remover will get rid of any files in patients folder ## 
	echo 'the remover.sh will go here'
	echo $i 'is a patient'
	#####################################
	## this will go through each time point for each patient ## 
	for a in {1..6}
	do
		echo $i $a 'json-to-fasta.sh'
		## this file should be open to append to, but I am curious to see how it 
		## does with multiple input files..## 
		echo $i $a 'get-v-genes.sh'		
	done;
	echo 'cat all 6 time points'
	echo 'make tsv'
	echo 'make ggplot'
	##############################################
	## this will align the gene and build the tree 
	for j in {1..7}
	do
		echo $i $j 'align gene'
		echo $i $j 'build tree with this aligned gene'
	done;
done;	
#########################################
	