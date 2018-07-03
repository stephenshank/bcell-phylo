
library("Glimma")
library("argparse")

#suppressPackageStartupMessages(library("argparse"))
args <- commandArgs
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args [2]

#in_dir <- args[1]

library(ggplot2)

d <- read.csv2(file, sep = '\t')

#d

ggplot(data = d) +
    labs(title = "77612", subtitle = "Clone size V-gene VS Visits, Sick at all time points.", x = 'Visits', y = 'Clone size of V-gene (log)') +
    geom_point(mapping = aes(Time,log(Size), color = Time, alpha =.03)) +
    facet_wrap(~V_gene) +
    theme(axis.text.x = element_text(face="bold", color="black", 
                           size=8, angle=45)) 
o_dir <- paste0(output,'/','V_usage.png')
ggsave(o_dir, width = 13, height =8)



