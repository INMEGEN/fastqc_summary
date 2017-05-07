library("ggplot2")
library("reshape2")
library("cowplot")

read_data <- function(file_name="summary.txt"){
	fastq <- read.table(file_name, header=FALSE)
	names(fastq)[1:8]<-c("Reads", "FileName", "Sequencer", "Run", "Flowcel", "Pair", "N", "Bad")
	fastq$filename <- gsub(fastq$FileName, pattern=".fq.gz", replacement="")
	fastq$filename <- gsub(fastq$FileName, pattern="data/", replacement="")
	fastq$Run <- factor(fastq$Run)
	design <- do.call(rbind, strsplit(gsub(fastq$FileName,
		pattern=".fq.gz", replacement=""), split="_"))
	colnames(design)<-c("Subject", "NA1", "Flowcell", "Lane", "PairEnd")
	fastq <- cbind(fastq, design)
	fastq$Subject <- factor(as.character(fastq$Subject), levels = sample(levels(fastq$Subject)))
	fastq$Flowcell <- factor(as.character(fastq$Flowcell), levels = sample(levels(fastq$Flowcell)))
	fastq$Lane <- as.integer(gsub(fastq$Lane, pattern="L", replacement=""))
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
	plot_design <- ggplot(data=fastq, aes(x=Lane, y=Subject, fill=Flowcell)) +
		geom_tile() + facet_grid(. ~ Sequencer) +
		theme(panel.grid.major = element_line("black", size = 0.1))
	return(plot_design)
}

plot_reads_per_subject <- function(fastq){
	total_reads <- tapply(fastq$Reads, INDEX=fastq$Subject, sum)
	total_reads <- sort(total_reads)
	subjects <- names(total_reads)
	fastq$Subject <- factor(as.character(fastq$Subject), levels=subjects)
	plot_reads <- ggplot(data=fastq, aes(x=Subject, y=log(Reads))) + 
		geom_boxplot() +
		geom_line(
			data=data.frame(
				 Subject=subjects, 
				Total_Reads=log(total_reads)
			), 
			aes(x=Subject, y=Total_Reads, group=1)
		)
	return(plot_reads)
}

fastq <- read_data()
fastq <- check_integrity(fastq)
#plot_sequencing(fastq)
plot_reads_per_subject(fastq)

#How many Flowcell were used by each subject?
table(rowSums(with(unique(fastq[, c("Subject", "Flowcell")]), table(Subject, Flowcell))))
# 1  2  3 
#63 19  4
#Are there sequencer effect? 
table(rowSums(with(unique(fastq[, c("Subject", "Sequencer")]), table(Subject, Sequencer))))
# 1  2  3 
#65 19  2 
