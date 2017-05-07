library("ggplot2")
library("reshape2")
library("cowplot")

read_data <- function(file_name="summary.txt"){
	fastq <- read.table(file_name, header=TRUE)
	names(fastq)[1]<-"reads"
	fastq$filename <- gsub(fastq$filename, pattern=".fq.gz", replacement="")
	fastq$filename <- gsub(fastq$filename, pattern="data/", replacement="")
	design <- do.call(rbind, strsplit(gsub(fastq$filename,
		pattern=".fq.gz", replacement=""), split="_"))
	colnames(design)<-c("Subject", "NA1", "Flowcell", "Lane", "PairEnd")
	fastq <- cbind(fastq, design)
	names(fastq)[1:8]<-c("Reads", "FileName", "Sequencer", "Run", "Flowcel", "Pair", "N", "Bad")
	return(fastq)
}

check_integrity <- function(fastq){
	stopifnot(all(fastq$Flowcel == fastq$Flowcell)); fastq$Flowcel <- NULL
	stopifnot(all(fastq$Pair == fastq$PairEnd)); fastq$PairEnd <- NULL
	stopifnot(all(fastq$Bad == 0)); fastq$Bad <- NULL
	stopifnot(all(fastq$N == "N")); fastq$N <- NULL
	return(fastq)
}

fastq <- read_data()
fastq <- check_integrity(fastq)

plot_design <- ggplot(data=fastq, aes(x=Lane, y=PairEnd, fill=Subject)) +
	geom_tile() + facet_grid(Sequencer ~ Run)
plot_design



