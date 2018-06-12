## this pipeline will run the modules one by one if need be. this was just cleaner than having one large bash file.
## eventually this will be a snakemake file or something similar

#lets put something here that will jump into these folders and erase files already there

#bash modules/json_to_fasta.sh

#bash modules/cat_unaligned.sh

#bash modules/get_v-genes.sh

#bash modules/aligner.sh
##############################################
## input = unaligned_jsons, unaligned_v-genes
## output = aligned_jsons, aligned_v-genes 
############################################

#bash modules/tsv_maker.sh
#############################
## this doesn't work right now
#############################

#bash modules/tree-builder.sh

#bash modules/ggplotter.sh

## need to make an 'erase' script, although the git ignore file might be sufficient


