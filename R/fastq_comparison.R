source("per_base_sequence_quality.R")
library("ggplot2")
library("cowplot")

raw <- quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")
bgi <- quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")
inmegen <- quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")

comparison <- plot_grid(
	plotlist=list(NULL, NULL, NULL,raw, bgi, inmegen),
	labels=c("Raw", "BGI", "INMEGEN"),
	ncol=3,
	nrow=2,
	rel_heights=c(0.05, 0.95)
)

comparison

