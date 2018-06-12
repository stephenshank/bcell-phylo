## this pipeline will run the modules one by one if need be. this was just cleaner than having one large bash file.
## eventually this will be a snakemake file or something similar

bash modules/remover.sh

bash modules/json_to_fasta.sh
##################################
## input = jsons from tar ball 
## output = unaligned fasta files 
#################################

bash modules/cat_unaligned.sh
#################################################################
## input = all unaligned json timepoint fasta files
## output = unaligned fasta files with all time points together 
################################################################

bash modules/get_v-genes.sh
######################################################
## input = full unaligned time-point json fasta file 
## output = each unique v gene unaligned fasta files
#####################################################

bash modules/aligner.sh
##############################################
## input = unaligned_jsons, unaligned_v-genes
## output = aligned_jsons time points, aligned_v-genes across all time points**
############################################

bash modules/tsv_maker.sh
#############################
## input = full aligned json 
## output = tsv file 
#############################

bash modules/tree-builder.sh
#############################################
## input = all aligned all jsons and v-genes
## output = trees in newish format 
############################################

bash modules/ggplotter.sh
##################################
## input = tsv file made earlier
## output = pretty png files
#################################
