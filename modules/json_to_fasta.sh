for i in {7..12}
do
  python python/vfull_json_to_fasta.py \
	--file data/sandbox/input_json/$i\_clone.json \
 	--size 30 \
	--output data/sandbox/output/unaligned_jsons/;
done; 
