rm(list=ls())
load("/home/cfresno/Dropbox/inmegen/100g/fq/R/fastqc_data.RData")
load(file="/home/cfresno/Dropbox/inmegen/100g/fq/R/Reads.RData")
rm(list=ls()[!ls()%in%c("raw_data", "bgi_data", "inmegen_data")])

source("per_base_sequence_quality.R")
source("per_sequence_quality_scores.R")
library("ggplot2")
library("reshape2")
library("cowplot")



##Read comparison
raw_psqs <- readPSQS("data/raw", group="Raw")
bgi_psqs <- readPSQS("data/clean", group="BGI")
inmegen_psqs <- readPSQS("data/trimmomatic", group="INMEGEN")
#save.image(file="results/Reads.RData", compress="xz")

#mean_quality <- plot_mean_quality_density(data)

#save.image(file="results/Comparison.RData", compress="xz")

