## this is the full list, we will just use the first two to try ##
#ids=('28729' '48689' '67029' '77612' '78202' '93954' '99361' '99682' 'GJS')
ids=('28729' '48689')

for i in ${ids[@]}
do 	
	rm -r data/out/$i
	## this remover will get rid of any files in patients folder ## 
	echo 'removing populated folders for' $i
	rm data/out/$i/$i*.tsv
	rm data/out/$i/*.fasta
	rm data/out/$i/$i*.new
	rm data/out/pretty_pictures/*.png
	
done;