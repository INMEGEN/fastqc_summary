options(width=120)
library("ggplot2")
library("reshape2")
library("cowplot")
library("splines")
library("plyr")
################################################################################
readPSQS <- function(data_dir, group="Raw"){
	files <- list.files(path=data_dir, pattern="Per_sequence_quality_scores", full.names=TRUE)
	data <- lapply(files, read.table, sep="\t", comment.char ="*", header=TRUE)
	data <- lapply(data,function(x){
        aux <- t(x)
        colnames(aux)<-aux[1,]
        return(aux[-1, , drop=FALSE])
	})
    data <- do.call(rbind.fill,lapply(data, as.data.frame))
    data <- colSums(data, na.rm=TRUE)
    return(data.frame(Quality=names(data), Count=data, Group=group))
}

plot_mean_quality_density<-function(data){
	data$Quality <- as.integer(as.character(data$Quality))
	x_breaks <- 2:42
	plot_obj <- ggplot(
        	data=data, 
            	aes(x=Quality, y=Count/10^6, group=Group, color=Group)) +
        #geom_line()+ 
        geom_smooth(se=FALSE)+
        xlab("Mean Sequence Quality (Phred Score)") +
        ylab("Number of reads [x 10^6]")+
        theme(legend.position="bottom")+
	scale_x_continuous(breaks=x_breaks, labels=x_breaks)
    return(plot_obj)
}
