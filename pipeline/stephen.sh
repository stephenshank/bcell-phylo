rm data/no_replicates/*.fasta
rm data/no_replicates/*.new
rm data/replicates/*.fasta
rm data/replicates/*.new

bash pipeline/replicates.sh
bash pipeline/no_replicates.sh
