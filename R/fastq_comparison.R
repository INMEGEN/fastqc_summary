rm(list=ls())
load("/home/cfresno/Dropbox/inmegen/100g/fq/R/fastqc_data.RData")
rm(list=ls()[!ls()%in%c("raw_data", "bgi_data", "inmegen_data")])

source("per_base_sequence_quality.R")
source("per_sequence_quality_scores.R")
library("ggplot2")
library("reshape2")
library("cowplot")

#raw_data <- read_data("data/raw")
#bgi_data <- read_data("data/clean")
#inmegen_data <- read_data("data/trimmomatic")

raw <- quality_plot(raw_data)
bgi <- quality_plot(bgi_data)
inmegen <- quality_plot(inmegen_data)

##Full data fastqc plots
comparison <- plot_grid(
	plotlist=list(raw, bgi, inmegen),
	labels=c("Raw", "BGI", "INMEGEN"),
	ncol=3,
	nrow=1,
	rel_heights=c(0.05, 0.95),
	hjust=-1
)
comparison
##FIXME to save the file
#ggsave(comparison, file="comparison.png", dpi=300, device="cairo")
#png(file="comparison.png", res=300, width=1024, height=720, type="cairo")
#comparison
#dev.off()
#Save with screenshoot, transparency not supported

## Summary fastqc plots
compare <- plot_grid(
	plotlist=list(
		fastqc_plot(fastq_summary(raw_data)),
		fastqc_plot(fastq_summary(bgi_data)),
		fastqc_plot(fastq_summary(inmegen_data))
	),
	labels=c("RAW", "BGI", "INMEGEN"),
	ncol=3,
	nrow=1,
	hjust=-1
)
compare

