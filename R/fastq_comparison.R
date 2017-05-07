rm(list=ls())
#load("/home/cfresno/ssh/castillo/tmp/Comparison.RData")
rm(list=ls()[!ls()%in%c("raw_data", "bgi_data", "inmegen_data")])

source("./R/per_base_sequence_quality.R")
source("./R/per_sequence_quality_scores.R")
library("ggplot2")
library("reshape2")
library("cowplot")

#data_dir_example <- "/home/cfresno/ssh/castillo/tmp/pbsq"

raw_data <- read_data("data/raw")
bgi_data <- read_data("data/clean")
inmegen_data <- read_data("data/trimmomatic")

raw <- quality_plot(raw_data)
bgi <- quality_plot(bgi_data)
inmegen <- quality_plot(inmegen_data)

##Full data fastqc plots
comparison <- plot_grid(
	plotlist=list(NULL, NULL, NULL,raw, bgi, inmegen),
	labels=c("Raw", "BGI", "INMEGEN"),
	ncol=3,
	nrow=2,
	rel_heights=c(0.05, 0.95)
)
comparison

complete<-rbind(
    cbind(raw_data, Group="RAW"),
    cbind(bgi_data, Group="BGI"),
    cbind(inmegen_data, Group="INMEGEN"))

complete<-cbind(inmegen_data, Group="INMEGEN")

fastqc_complete <- function(complete){
    ##ggplot compatible data
    complete_melt <- melt(complete[, c("Base", "Mean", "Q2", "Q1", "Q3", "P10", "P90", "x", "Group")],
        id=c("Base", "x", "Group"), variable.name = "Measurement", value.name = "Quality")
    complete_melt$Region <- factor(sapply(as.character(complete_melt$Measurement), function(x){
        switch(x, 
            "Mean"= "Mean-Median", "Q2"="Mean-Median",
            "Q1"="Q1-Q3", "Q3"="Q1-Q3", 
            "P10"="P10-P90", "P90"="P10-P90")
    }))
    complete_melt$Region<-factor(as.character(complete_melt$Region),
        levels=c("P10-P90", "Mean-Median", "Q1-Q3"))
    #head(complete_melt)
    ##Polygon regions
    region_coord <- lapply(levels(complete$Group), function(group){
        datum <- subset(complete, Group==group)
        rbind(
            cbind(polygon_from_region(datum$P10, datum$P90, datum$x), Group=group, Region="P10-P90"),
            cbind(polygon_from_region(datum$Q1, datum$Q3, datum$x), Group=group, Region="Q1-Q3"),
            cbind(polygon_from_region(datum$Mean, datum$Q2, datum$x), Group=group, Region="Mean-Median")
            )
    })
    
    ##Colors
    region_colors <- c("P10-P90"="#eded44", "Q1-Q3"="#01cbf3", "Mean-Median"="#43640b")
    mean_color<-"#fdd65d"
    median_color<-"#f6a801"

    ##Axis variables
    x_breaks <- c(1:9, seq(from=11, to=38, by=3))
    x_labels <- levels(complete_melt$Base)[x_breaks]
    y_breaks <- seq(from=0, to=42, by=2)
    
    ##The plot
    p <- ggplot(data=complete_melt, aes(x=x, y=Quality, color=Region, Group=Measurement)) + 
        facet_grid(. ~ Group) + 
        geom_point(alpha=0.3) + geom_jitter(alpha=0.3)+
        ##Regions
        geom_polygon(data=region_coord,
        aes(x=x,y=y),
        fill=region_color,
        alpha=0.1
    )
        ##Quality_thresholds
        geom_hline(aes(yintercept=20), linetype="dashed", color="red", alpha=0.5) +
        geom_hline(aes(yintercept=28), linetype="dashed", color="green", alpha=0.5)+
        ##Axis
        scale_x_continuous(name="Position in read (bp)", breaks=x_breaks, labels=x_labels) +
        scale_y_continuous(name=, "Quality (Phred Score)", breaks=y_breaks, 
            labels=as.character(y_breaks), limits=c(0, 42))+
        ##Color scale
        scale_color_manual(values=region_colors, drop=FALSE)+    
        ##Theme
        theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1))
      
    p   
    
    
}


## Summary fastqc plots
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

#data_dir<-"/home/cfresno/ssh/castillo/tmp/summary/data"
#data <- readPSQS(data_dir, group="Raw")
#mean_quality <- plot_mean_quality_density(data)

raw_psqs <- readPSQS("data/raw", group="Raw")
bgi_psqs <- readPSQS("data/clean", group="BGI")
inmegen_psqs <- readPSQS("data/trimmomatic", group="INMEGEN")

mean_quality <- plot_mean_quality_density(data)

#save.image(file="results/Comparison.RData", compress="xz")
