source("per_base_sequence_quality.R")
library("ggplot2")
library("reshape2")
library("cowplot")

raw <- quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")
bgi <- quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")
inmegen <- quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")

##Full data fastqc plots
comparison <- plot_grid(
	plotlist=list(NULL, NULL, NULL,raw, bgi, inmegen),
	labels=c("Raw", "BGI", "INMEGEN"),
	ncol=3,
	nrow=2,
	rel_heights=c(0.05, 0.95)
)

comparison

## Summary fastqc plots
raw_data <- read_data("/home/cfresno/ssh/castillo/tmp/pbsq")

compare <- plot_grid(
	plotlist=list(NULL, NULL, 
        fastqc_plot(fastq_summary(raw_data)),
        fastqc_summary_boxplot(raw_data)),
	labels=c("Summary", "Using boxplots"),
	ncol=2,
	nrow=2,
	rel_heights=c(0.05, 0.95),
	hjust=-1
)
compare





