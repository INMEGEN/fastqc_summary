library("ggplot2")
library("reshape2")
library("cowplot")

#function(name)
#	names(fastq)<-sapply(names(fastq), function(x){
#		paste(toupper(substr(x,1,1)),
#			substr(x, 2, nchar(x)), sep="")
#	})
#}

read_data <- function(file_name="experiment.stats"){
	fastq <- read.table(file_name, header=TRUE)
	fastq$Run <- factor(fastq$Run)
	design <- design_from_filename(fastq$Filename)
	colnames(design)<-c("Subject", "NA1", "Flowcel", "Lanes", "PairEnd")
	fastq <- cbind(fastq, design)
	fastq$Subject <- factor(as.character(fastq$Subject), levels = sample(levels(fastq$Subject)))
	fastq$Flowcel <- factor(as.character(fastq$Flowcel), levels = sample(levels(fastq$Flowcel)))
	fastq$Lanes <- gsub(fastq$Lanes, pattern="L", replacement="")
	return(fastq)
}

design_from_filename <- function(filename) {
	basename <- gsub(filename, pattern=".fq.gz", replacement="")
	design <- strsplit(basename, split="_")
	design <- do.call(rbind, design)
	return(design)
}

check_integrity <- function(fastq){
	stopifnot(all(fastq$Flowcel == fastq$Flowcell)); fastq$Flowcel <- NULL
	stopifnot(all(fastq$Pair == fastq$PairEnd)); fastq$PairEnd <- NULL
	stopifnot(all(fastq$Lane == fastq$Lanes)); fastq$Lanes <- NULL
	return(fastq)
}

##FIXME change the grid to border each subject
plot_sequencing <- function(fastq){
	plot_design <- ggplot(data=fastq, aes(x=Lane, y=Subject, fill=Flowcell)) +		geom_tile() + facet_grid(. ~ Sequencer) +
		theme(panel.grid.major = element_line("black", size = 0.1))+
		guides(fill=guide_legend(ncol=1))
	return(plot_design)
}

plot_reads_per_subject <- function(fastq){
	total_reads <- tapply(fastq$Reads, INDEX=fastq$Subject, sum)
	total_reads <- sort(total_reads)
	subjects <- names(total_reads)
	fastq$Subject <- factor(as.character(fastq$Subject), levels=subjects)
	plot_reads <- ggplot(data=fastq, aes(x=Subject, y=Reads/10^6)) + 
		geom_boxplot() +
		scale_x_discrete("Subject", labels=1:length(subjects))+ 
		ylim(0, 120)+
		ylab("Reads [x 10^6]") +
		geom_hline(aes(yintercept=mean(fastq$Reads)/10^6, color="Mean"), linetype="dashed") +
		scale_colour_manual("", values=c("Mean"="blue"))
	return(plot_reads)
}

fastq <- read_data("../data/raw.info")
fastq <- check_integrity(fastq)
sequencing <- plot_sequencing(fastq)
sequencing
#Screenshot saving plot

#Reads level
raw <- check_integrity(read_data("../data/raw.info"))
bgi <- check_integrity(read_data("../data/clean.info"))
inmegen_paired <- check_integrity(read_data("../data/trimmomatic.paired.info"))
inmegen_unpaired <- check_integrity(read_data("../data/trimmomatic.unpaired.info"))

inmegen <- inmegen_paired
inmegen$Reads <- inmegen_paired$Reads + inmegen_unpaired$Reads

raw_plot <- plot_reads_per_subject(raw)
bgi_plot <- plot_reads_per_subject(bgi)
inmegen_plot <- plot_reads_per_subject(inmegen)

reads <- plot_grid(
	plotlist=list(raw_plot, bgi_plot, inmegen_plot),
	ncol=3,
	nrow=1,
	labels=c("RAW", "BGI", "INMEGEN"),
	hjust=-1
)
reads

#How many Flowcell were used by each subject?
#table(rowSums(with(unique(fastq[, c("Subject", "Flowcell")]), table(Subject, Flowcell))))
# 1  2  3 
#63 19  4
#Are there sequencer effect? 
#table(rowSums(with(unique(fastq[, c("Subject", "Sequencer")]), table(Subject, Sequencer))))
# 1  2  3 
#65 19  2 
