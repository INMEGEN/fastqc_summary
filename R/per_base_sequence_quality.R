################################################################################
##Get the summary of multiple fastqc_data.txt.* results:
## -Per_base_sequence_quality:
##Author: Cristobal Fresno & Joshua Haase
##Date: 2017/05/01
################################################################################
#options(width=80)
options(width=140)
library("ggplot2")
library("reshape2")
library("cowplot")
library("splines")
################################################################################
##-Per_base_sequence_quality
################################################################################
readPBSQ <- function(filePBSQ){
    pbsq <- lapply(filePBSQ, read.table, sep="\t", comment.char ="*",
        header=TRUE)
    pbsq <- do.call(rbind,lapply(1:length(filePBSQ), function(x){
        pbsq[[x]]$File <- filePBSQ[x]
        return(pbsq[[x]])
    }))
    names(pbsq) <- c("Base", "Mean", "Q2", "Q1", "Q3", "P10", "P90", "File")

    ##Change the name to include additional data.
    aux <- as.data.frame(do.call(rbind, strsplit(gsub(pbsq$File,
        pattern="_fastqc_data.txt.Per_base_sequence_quality", replacement=""),
        split="_")))
    names(aux) <- c("Subject", "NA1", "Flowcell", "Lane", "PairEnd")
    pbsq$File <- NULL
    pbsq <- cbind(pbsq, aux)
    ##Order the Base levels
    pbsq$Base <- factor(as.character(pbsq$Base),
            levels=levels(pbsq$Base)[
                order(as.integer(unlist(lapply(strsplit(
                    levels(pbsq$Base), split="-"), function(x)x[1]))))])
    ##Clean the Subject Path
    aux <- do.call(rbind, strsplit(as.character(pbsq$Subject), split="/"))
    pbsq$Subject <- aux[,ncol(aux)]
    return(pbsq)
}

#Transparency work around
# http://tinyheero.github.io/2015/09/15/semi-transparency-r.html
#in ~./Rprofile
#setHook(packageEvent("grDevices", "onLoad"),
#function(...) grDevices::X11.options(type='cairo'))
#options(device='x11')

smooth_polygon_edge <- function(coordinates) {
    smooth_model <- loess(y ~ x, data=coordinates)
    smooth_data <- predict(smooth_model)
    return(smooth_data)
}

polygon_from_region <- function(inferior_limit, superior_limit, x){
    xMax <- max(x)
    stopifnot(length(inferior_limit) == length(superior_limit))
    min_down <- tapply(inferior_limit, INDEX=x, min)
    max_top <- tapply(superior_limit, INDEX=x, max)
    inf_line <- smooth_polygon_edge(data.frame(x=1:xMax, y=min_down))
    sup_line <- smooth_polygon_edge(data.frame(x=1:xMax, y=max_top))
    return(
        data.frame(
            x=c(1:xMax, xMax:1),
            y=c(inf_line, sup_line[xMax:1])
        )
    )
}

plot_region <- function(inferior_limit, superior_limit, x, region_color, dots_color, plot=ggplot()) {
    data <- data.frame(inferior_limit, superior_limit, x)
    return(
    plot+
    geom_polygon(
        data=polygon_from_region(inferior_limit, superior_limit, x),
        aes(x=x,y=y),
        fill=region_color,
        alpha=0.1
    )+
    geom_jitter(data=data, aes(x=x, y=superior_limit), color=dots_color, alpha=0.3)+
    geom_jitter(data=data, aes(x=x, y=inferior_limit), color=dots_color, alpha=0.3))
}

plot_quality_limits <- function(plot) {
    return(plot +
        geom_hline(aes(yintercept=20), linetype="dashed", color="red", alpha=0.5) +
        geom_hline(aes(yintercept=28), linetype="dashed", color="green", alpha=0.5))
}

plot_trend_line <- function(coordinates, color, plot=ggplot()) {
    return(plot +
        geom_smooth(data=coordinates, aes(x=x, y=y), color=color, se=FALSE))
}

read_data <- function(data_dir) {
	filePBSQ <- list.files(path=data_dir, pattern="Per_base_sequence_quality", full.names=TRUE)
	# per base sequence quality
	pbsq <- readPBSQ(filePBSQ)
	pbsq$x <- as.integer(unclass(pbsq$Base))
	return(pbsq)
}

plot_labels <-function(plot, x_label, y_label) {
	plot <- plot + xlab(x_label)
	plot <- plot + ylab(y_label)
	return(plot)
}

plot_legend <-function(plot) {
	out_color="#eded44"
	middle_color="#01cbf3"
	in_color="#43640b"
	mean_color="#fdd65d"
        median_color="#f6a801"

	labels <- data.frame(Regions=c(out_color, middle_color, in_color), labels=c("P10-P90", "Q1-Q3", "Mean-Median"))
	values <- labels$Regions; names(values)=labels$labels
	glegend <- ggplot(data=labels, aes(x=1, y=1, fill=Regions, label=labels))+
		geom_tile() +
		scale_fill_manual(name="Regions", values=values)
	glegend

    x_breaks <- c(1:9, seq(from=11, to=38, by=3))
	x_labels <- levels(data$Base)[x_breaks]
	plot <- plot + scale_x_continuous(breaks=x_breaks, labels=x_labels)
	y_breaks <- seq(from=0, to=42, by=2)
	plot <- plot + scale_y_continuous(breaks=y_breaks, labels=as.character(y_breaks))
	return(plot)
}

plot_axis <-function(plot_obj, data) {
	x_breaks <- c(1:9, seq(from=11, to=38, by=3))
	x_labels <- levels(data$Base)[x_breaks]
	plot_obj <- plot_obj + scale_x_continuous(breaks=x_breaks, labels=x_labels)
	y_breaks <- seq(from=0, to=42, by=2)
	plot_obj <- plot_obj + scale_y_continuous(breaks=y_breaks, labels=as.character(y_breaks), limits=c(0,42))
	return(plot_obj)
}

quality_plot <- function(data, out_color="#eded44", middle_color="#01cbf3", in_color="#43640b", mean_color="#fdd65d", median_color="#f6a801") {
	quality_p <- plot_region(data$P10, data$P90, data$x, out_color, out_color)
	quality_p <- plot_region(data$Q1, data$Q3, data$x, middle_color, middle_color, quality_p)
	quality_p <- plot_region(data$Mean, data$Q2, data$x, in_color, in_color, quality_p)
	quality_p <- plot_quality_limits(quality_p)
	quality_p <- plot_trend_line(coordinates=data.frame(x=data$x, y=data$Q2), median_color, quality_p)
	quality_p <- plot_trend_line(coordinates=data.frame(x=data$x, y=data$Mean), mean_color, quality_p)
	quality_p <- plot_labels(quality_p, "Position in read (bp)", "Quality (10 ⨉ -log(pe))")
	quality_p <- plot_axis(quality_p, data)
	return(quality_p)
}

##Single per base sequence quality plot#########################################
fastqc_plot <-  function(datum){
    datum$Base <- factor(as.character(datum$Base),
        levels=levels(datum$Base)[
            order(as.integer(unlist(lapply(strsplit(
                levels(datum$Base), split="-"), function(x)x[1]))))])
    datum$x <- 1:nrow(datum)
    xMax <- length(levels(datum$Base))
    background <- data.frame(
                    x=rep(c(0,xMax,xMax,0),3),
                    y=c(c(0,0,20,20),c(20,20,28,28),c(28,28,42,42)),
                    Quality=c(rep("Bad", 4),rep("Intermediate", 4),rep("Good", 4)))
    sp <- ggplot()+
        geom_polygon(data=background, aes(x=x, y=y, group=Quality, fill=Quality))+
        geom_boxplot(data=datum, aes(x=x, ymin=P10, lower=Q1, middle=Q2,
            upper=Q3, ymax=P90, group=Base), stat = "identity", fill="yellow")+
        geom_line(data=datum, aes(x=x, y=Mean, group=1), color="blue")+
        theme(legend.position="bottom")
    sp <- plot_axis(sp, datum)
    sp <- plot_labels(sp, "Position in read (bp)", "Quality (10 ⨉ -log(pe))")
    return(sp)
}

fastqc_summary_boxplot <-function(data){
    fastqc_melt <- melt(data[, c(1:7, ncol(data))], id=c("Base", "x"))
    mean_points <- tapply(fastqc_melt$value, INDEX=fastqc_melt$Base, mean)
    data_mean <- data.frame(x=names(mean_points), Mean=mean_points)
    x_breaks <- c(1:9, seq(from=11, to=38, by=3))
    x_labels <- levels(data$Base)
    x_labels[-x_breaks] <- ""
    y_breaks <- seq(from=0, to=42, by=2)
    xMax <- max(fastqc_melt$x)
    background <- data.frame(
                    x=rep(c(0,xMax,xMax,0),3),
                    y=c(c(0,0,20,20),c(20,20,28,28),c(28,28,42,42)),
                    Quality=c(rep("Bad", 4),rep("Intermediate", 4),rep("Good", 4)))
    p <- ggplot(fastqc_melt, aes(x=Base, y=value)) + 
        geom_polygon(data=background, aes(x=x, y=y, group=Quality, fill=Quality))+
        geom_boxplot(fill="yellow") +
        geom_line(data=data_mean, aes(x=x, y=Mean, group=1), color="blue")+
        theme(legend.position="bottom")
    p <- plot_labels(p, "Position in read (bp)", "Quality (10 ⨉ -log(pe))")
    p <- p + scale_x_discrete(labels=x_labels)
    p <- p + ylim(0, 42) + scale_y_continuous(breaks=y_breaks, labels=as.character(y_breaks))
    return(p)
}

fastq_summary <-function(data){
    return(data.frame(
        Base = as.factor(levels(data$Base)),
        P90 = with(data, tapply(data$P90, INDEX=data$Base, max)),
        Q3 = with(data, tapply(data$Q3, INDEX=data$Base, median)),
        Q2 = with(data, tapply(data$Q2, INDEX=data$Base, min)),
        Mean = with(data, tapply(data$Mean, INDEX=data$Base, min)),
        Q1 = with(data, tapply(data$Q1, INDEX=data$Base, min)),
        P10 = with(data, tapply(data$P10, INDEX=data$Base, min))))
}

test <- function() {
    data <- read_data("/home/cfresno/ssh/castillo/tmp/pbsq")
	datum <- subset(data, Subject=="SM-3MG3L" & Lane=="L1" & PairEnd=="1")
	fastqc_plot(datum)
	quality_plot("/home/cfresno/ssh/castillo/tmp/pbsq")
	fastqc_plot(fastq_summary(data))
}
