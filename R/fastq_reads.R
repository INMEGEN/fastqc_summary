rm(list=ls())
load("/home/cfresno/Dropbox/inmegen/100g/fq/R/Reads.RData")
rm(list=ls()[!ls()%in%c("raw_psqs", "bgi_psqs", "inmegen_psqs")])

source("per_base_sequence_quality.R")
source("per_sequence_quality_scores.R")
library("ggplot2")
library("reshape2")
library("cowplot")



##Read comparison
# raw_psqs <- readPSQS("data/raw", group="Raw")
# bgi_psqs <- readPSQS("data/clean", group="BGI")
# inmegen_psqs <- readPSQS("data/trimmomatic", group="INMEGEN")
# save.image(file="results/Reads.RData", compress="xz")

plot_mean_quality_density(rbind(raw_psqs, bgi_psqs, inmegen_psqs))

