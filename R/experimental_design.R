library("ggplot2")
library("reshape2")
library("cowplot")

read_data <- function(file_name="summary.txt"){
	fastq <- read.table(file_name, header=FALSE)
	names(fastq)[1:8]<-c("Reads", "FileName", "Sequencer", "Run", "Flowcel", "Pair", "N", "Bad")
	fastq$filename <- gsub(fastq$FileName, pattern=".fq.gz", replacement="")
	fastq$filename <- gsub(fastq$FileName, pattern="data/", replacement="")
	design <- do.call(rbind, strsplit(gsub(fastq$FileName,
		pattern=".fq.gz", replacement=""), split="_"))
	colnames(design)<-c("Subject", "NA1", "Flowcell", "Lane", "PairEnd")
	fastq <- cbind(fastq, design)
	return(fastq)
}

check_integrity <- function(fastq){
	stopifnot(all(fastq$Flowcel == fastq$Flowcell)); fastq$Flowcel <- NULL
	stopifnot(all(fastq$Pair == fastq$PairEnd)); fastq$PairEnd <- NULL
	stopifnot(all(fastq$Bad == 0)); fastq$Bad <- NULL
	stopifnot(all(fastq$N == "N")); fastq$N <- NULL
	return(fastq)
}

plot_sequencing <- function(fastq){
	plot_design <- ggplot(data=fastq, aes(x=Lane, y=Pair, fill=Subject)) +
		geom_tile() + facet_grid(Run ~ Sequencer)
	plot_design
	return(plot_design)
}

fastq <- read_data()
fastq <- check_integrity(fastq)
plot_sequencing(fastq)



