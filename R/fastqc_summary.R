################################################################################
##Get the summary of multiple fastqc_data.txt.* results:
## -Basic_Statistics: File, Type, Encoding, Sequences, PoorQC, Length, %GC.
## -Sequence_Length_Distribution: Are there different distributions pre subject?
## -Overrepresented_sequences: Are there any adapters?
## -Per_base_sequence_quality:
## -Adapter_Content
## -Kmer_Content
## -Per_base_N_content
## -Per_base_sequence_content
## -Per_sequence_GC_content
## -Per_sequence_quality_scores
## -Per_tile_sequence_quality
## -Sequence_Duplication_Levels
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
##-Basic_Statistics
##  File, Type, Encoding, Sequences, PoorQC, Length, %GC
################################################################################
fileBS<-list.files(path=".", pattern="Basic_Statistics")
bs<-lapply(fileBS, read.table, sep="\t")
bs<-as.data.frame(do.call(rbind,lapply(bs, function(x){t(x)[-1,]})))
names(bs)<-c("File", "Type", "Encoding", "Sequences", "PoorQC", "Length", "%GC")
bs
#                                     File                    Type              Encoding Sequences PoorQC Length %GC
# 1 SM-3MG3L_DHG02288_H3VGKCCXX_L1_1.fq.gz Conventional base calls Sanger / Illumina 1.9  38251775      0    150  44
# 2 SM-3MG3L_DHG02288_H3VGKCCXX_L1_2.fq.gz Conventional base calls Sanger / Illumina 1.9  38251775      0    150  43
################################################################################
##-Sequence_Length_Distribution
##  Are there different distributions per subject?
################################################################################
fileSLD<-list.files(path=".", pattern="Sequence_Length_Distribution")
sld<-lapply(fileSLD, read.table, sep="\t", comment.char ="*", header=TRUE)
sld<-do.call(rbind,lapply(1:length(fileSLD), function(x){
    sld[[x]]$File<-fileSLD[x]
    return(sld[[x]])
}))
names(sld)<-gsub(pattern="X.", replacement="", names(sld))
sld
#   Length    Count                                                                          File
# 1    150 38251775 SM-3MG3L_DHG02288_H3VGKCCXX_L1_1_fastqc_data.txt.Sequence_Length_Distribution
# 2    150 38251775 SM-3MG3L_DHG02288_H3VGKCCXX_L1_2_fastqc_data.txt.Sequence_Length_Distribution

################################################################################
##-Overrepresented_sequences
##  Are there different distributions per subject?
################################################################################
fileOS<-list.files(path=".", pattern="Overrepresented_sequences")
os<-do.call(rbind,lapply(fileOS, read.table, sep="\t", comment.char ="*", 
    header=TRUE))
names(os)<-gsub(pattern="X.", replacement="", names(os))
totalOS<-cbind(unique(os[,c("Sequence", "Possible.Source")]),  
    Total=tapply(os$Count, INDEX=os$Possible.Source, sum))
totalOS
#                                             Sequence                                   Possible.Source  Total
# 1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGC         TruSeq Adapter, Index 11 (100% over 50bp) 231348
# 2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG Illumina Single End PCR Primer 1 (100% over 50bp) 249096    

################################################################################
##-Per_base_sequence_quality
################################################################################
readPBSQ<-function(filePBSQ){
    pbsq<-lapply(filePBSQ, read.table, sep="\t", comment.char ="*", 
        header=TRUE)
    pbsq<-do.call(rbind,lapply(1:length(filePBSQ), function(x){
        pbsq[[x]]$File<-filePBSQ[x]
        return(pbsq[[x]])
    }))
    names(pbsq)<-c("Base", "Mean", "Q2", "Q1", "Q3", "P10", "P90", "File")

    ##Change the name to include additional data.
    aux<-as.data.frame(do.call(rbind, strsplit(gsub(pbsq$File, 
        pattern="_fastqc_data.txt.Per_base_sequence_quality", replacement=""),
        split="_")))
    names(aux)<-c("Subject", "NA1", "Flowcell", "Lane", "PairEnd")
    pbsq$File<-NULL
    pbsq<-cbind(pbsq, aux)
    ##Order the Base levels
    pbsq$Base<-factor(as.character(pbsq$Base), 
            levels=levels(pbsq$Base)[
                order(as.integer(unlist(lapply(strsplit(
                    levels(pbsq$Base), split="-"), function(x)x[1]))))])
    return(pbsq)
}
filePBSQ<-list.files(path=".", pattern="Per_base_sequence_quality")
pbsq<-readPBSQ(filePBSQ)

##Single per base sequence quality plot#########################################
singlePerBaseQuality<- function(datum){
    datum$Base<-factor(as.character(datum$Base), 
        levels=levels(datum$Base)[
            order(as.integer(unlist(lapply(strsplit(
                levels(datum$Base), split="-"), function(x)x[1]))))])
    datum$x<-1:nrow(datum)    
    xMax<-length(levels(datum$Base))
    background<-data.frame(
                    x=rep(c(0,xMax,xMax,0),3),
                    y=c(c(0,0,20,20),c(20,20,28,28),c(28,28,42,42)),
                    Quality=c(rep("Bad", 4),rep("Intermediate", 4),rep("Good", 4)))
    sp<-ggplot()+
        geom_polygon(data=background, aes(x=x, y=y, group=Quality, fill=Quality))+
        geom_boxplot(data=datum, aes(x=x, ymin=P10, lower=Q1, middle=Q2, 
            upper=Q3, ymax=P90, group=Base), stat = "identity", fill="yellow")+
        geom_line(data=datum, aes(x=x, y=Mean, group=Subject), color="blue")+
        xlab("Position in read (bp)")+
        scale_x_continuous(breaks=datum$x, labels=datum$Base)+
        theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="bottom")       
    return(sp)
}
  
datum<-subset(pbsq, Subject=="SM-3MG3L" & Lane=="L1" & PairEnd=="1")
singlePerBaseQuality(datum)

##Smooth per base sequence quality plot#########################################
#Transparency work around 
# http://tinyheero.github.io/2015/09/15/semi-transparency-r.html
#in ~./Rprofile
#setHook(packageEvent("grDevices", "onLoad"),
#function(...) grDevices::X11.options(type='cairo'))
#options(device='x11')

data<-pbsq
#smoothPerBaseQuality<- function(data){
    data
    type<-"Mean"
    data$x<-as.integer(unclass(data$Base))
    xMax<-length(levels(data$Base))
    background<-data.frame(
                    x=rep(c(0,xMax,xMax,0),3),
                    y=c(c(0,0,20,20),c(20,20,28,28),c(28,28,42,42)),
                    Quality=c(rep("Bad", 4),rep("Intermediate", 4),rep("Good", 4)))
    tapply(data$P10, INDEX=data$x,sd)                
                
    sp<-ggplot()+
        geom_hline(aes(yintercept=20), linetype="dashed", color="red", alpha=0.5)+
        geom_hline(aes(yintercept=28), linetype="dashed", color="green", alpha=0.5)+
        geom_jitter(data=data, aes(x=x, y=P90), alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=P90), se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=P10), alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=P10), se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Q1), color="green", alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=Q1), color="green", se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Q3), color="green", alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=Q3), color="green", se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Q2), color="brown", alpha=0.03)+
        geom_smooth(data=data, aes(x=x, y=Q2), color="brown", se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Mean), color="purple", alpha=0.03)+
        geom_smooth(data=data, aes(x=x, y=Mean), color="purple", se=FALSE)+
        xlab("Position in read (bp)")+
        scale_x_continuous(breaks=data$x, labels=data$Base)+
        theme(axis.text.x=element_text(angle=45, hjust=1))
        ##FIXME: add color legend!!!
    sp     
        
    ##Region plot
    ##R1: min(P10) - max(P90)
    r1<-data.frame(
        x=c(1:xMax, xMax:1),
        y=c(tapply(data$P10, INDEX=data$x, min),
            tapply(data$P90, INDEX=data$x, max)[xMax:1]))

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

region1 <- polygon_from_region(data$P10, data$P90, data$x)
region2 <- polygon_from_region(data$Q1, data$Q3, data$x)
region3 <- polygon_from_region(data$Mean, data$Q2, data$x)

spr<-ggplot()+
    geom_polygon(data=region1, aes(x=x,y=y), fill="grey96", alpha=1)+
    geom_jitter(data=data, aes(x=x, y=P90), color="grey96", alpha=0.01)+
    geom_jitter(data=data, aes(x=x, y=P10), color="grey96", alpha=0.01)+
    geom_polygon(data=region2, aes(x=x,y=y), fill="grey72", alpha=1)+
    geom_jitter(data=data, aes(x=x, y=P90), color="grey72", alpha=0.01)+
    geom_jitter(data=data, aes(x=x, y=P10), color="grey72", alpha=0.01)+
    geom_polygon(data=region3, aes(x=x,y=y), fill="grey32", alpha=1)+
    geom_jitter(data=data, aes(x=x, y=P90), color="grey32", alpha=0.01)+
    geom_jitter(data=data, aes(x=x, y=P10), color="grey32", alpha=0.01)
    
spr

        geom_hline(aes(yintercept=20), linetype="dashed", color="red", alpha=0.5)+
        geom_hline(aes(yintercept=28), linetype="dashed", color="green", alpha=0.5)+
        geom_jitter(data=data, aes(x=x, y=P90), alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=P90), se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=P10), alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=P10), se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Q1), color="green", alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=Q1), color="green", se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Q3), color="green", alpha=0.03)+
        #geom_smooth(data=data, aes(x=x, y=Q3), color="green", se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Q2), color="brown", alpha=0.03)+
        geom_smooth(data=data, aes(x=x, y=Q2), color="brown", se=FALSE)+
        geom_jitter(data=data, aes(x=x, y=Mean), color="purple", alpha=0.03)+
        geom_smooth(data=data, aes(x=x, y=Mean), color="purple", se=FALSE)+
        xlab("Position in read (bp)")+
        scale_x_continuous(breaks=data$x, labels=data$Base)+
        theme(axis.text.x=element_text(angle=45, hjust=1))
        ##FIXME: add color legend!!!
    spr     
    
    head(data) 
    tapply(data$P10, INDEX=data$Base, median)
    
    
 











Histograma antes y despues a nivel del promedio de la lectura?
Per_base_sequence_quality
#Quality        Count
10      7.0
11      49.0
12      92.0
13      239.0
14      460.0


