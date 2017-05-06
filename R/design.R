################################################################################
##Get the experimental design out of the fq.gz files
##Author: Cristobal Fresno & Joshua Haase
##Date: 2017/04/30
################################################################################
options(width=80)
################################################################################
##Using file name description only
################################################################################
#Create the file name descritop fq.gz file 
#system("ls /100g/data/raw/*.fq.gz > raw.txt")
#Read the data
design<-read.table(file="raw.txt", sep=" ")
names(design)<-"File"
design$File<-as.character(design$File)
design<-cbind(design, 
    do.call(rbind, 
        strsplit(gsub(design$File, pattern=".fq.gz", replacement=""),
            split="_")))
names(design)<-c("FileName", "Subject", "NA1", "Flowcell", "Lane", "PairEnd")
dim(design)
# [1] 1158    6
head(design)
#                                 FileName  Subject      NA1  Flowcell Lane PE
# 1 SM-3MG3L_DHG02288_H3VGKCCXX_L1_1.fq.gz SM-3MG3L DHG02288 H3VGKCCXX   L1 1
# 2 SM-3MG3L_DHG02288_H3VGKCCXX_L1_2.fq.gz SM-3MG3L DHG02288 H3VGKCCXX   L1 2
# 3 SM-3MG3L_DHG02288_H3VGKCCXX_L2_1.fq.gz SM-3MG3L DHG02288 H3VGKCCXX   L2 1
# 4 SM-3MG3L_DHG02288_H3VGKCCXX_L2_2.fq.gz SM-3MG3L DHG02288 H3VGKCCXX   L2 2
# 5 SM-3MG3L_DHG02288_H3VGKCCXX_L3_1.fq.gz SM-3MG3L DHG02288 H3VGKCCXX   L3 1
# 6 SM-3MG3L_DHG02288_H3VGKCCXX_L3_2.fq.gz SM-3MG3L DHG02288 H3VGKCCXX   L3 2
#Total                                1158       95       95        25    8 2

##Exploring each variable
with(design,length(unique(as.character(Subject))))
# [1] 95
with(design,length(unique(as.character(NA1))))
# [1] 95
# Are the Subject and NA1 field coding the same data?   
# Yes. This is the extended subject name!!
contingecy<-with(design,table(Subject,NA1))
stopifnot(sum(diag(contingecy)) == sum(contingecy))

with(design,length(unique(as.character(Flowcell))))
# [1] 25
with(design,length(unique(as.character(Lane))))
# [1] 8
#Is there a preferential lane??
#No. That's OK.
with(design, table(Lane))
# Lane
#  L1  L2  L3  L4  L5  L6  L7  L8 
# 166 102 160 154 134 148 134 160 
summary(with(design, as.numeric(table(Lane))))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   102.0   134.0   151.0   144.8   160.0   166.0 

# Is the PairEnd per Lane design balanced?
# Yes. That's OK.
with(design,length(unique(as.character(PairEnd))))
# [1] 2
with(design,table(PairEnd, Lane))
# PairEnd L1 L2 L3 L4 L5 L6 L7 L8
#       1 83 51 80 77 67 74 67 80
#       2 83 51 80 77 67 74 67 80

##Crossing the different available variables
# Are the subjects confused with the Flowcells? 
# Yes!!! Some Flowcell only contains one subject. This design is not the 
# best one at all by far.
with(design,tapply(Subject, INDEX=Flowcell, function(x){length(unique(x))}))
# H2KH2CCXX H2KV5CCXX H2L7GCCXX H2L7JCCXX H2L7KCCXX H2L7WCCXX H2LC5CCXX H2LCVCCXX 
#         2         6         1         8         8         8         5         3 
# H3TKWCCXX H3TLVCCXX H3TVLCCXX H3V3GCCXX H3V3MCCXX H3VGCCCXX H3VGFCCXX H3VGKCCXX 
#         3         1         3        10         9         9         6         7 
# H3VJ2CCXX H3VJHCCXX H3VKYCCXX H3W3HCCXX H3WFHCCXX H3WG2CCXX H7GYTCCXX HF2GTCCXX 
#         5         2         8         3         1         1         9        26 
# HF5T2CCXX 
#         1

# How many different Flowcell each Subject has?
# Half of the subjects in a single Flowcell, is not a good idea.
with(design,table(tapply(Flowcell, INDEX=Subject, 
    function(x){length(unique(x))})))
#  1  2  3  4 
# 53 35  6  1 
with(design,tapply(Flowcell, INDEX=Subject, function(x){length(unique(x))}))
# SM-3MG3L SM-3MG3M SM-3MG3N SM-3MG3O SM-3MG3P SM-3MG3R SM-3MG3U SM-3MG3V 
#        1        1        1        2        2        1        2        2 
# SM-3MG3Y SM-3MG3Z SM-3MG45 SM-3MG46 SM-3MG47 SM-3MG48 SM-3MG49 SM-3MG4A 
#        1        2        2        2        2        2        1        1 
# SM-3MG4B SM-3MG4C SM-3MG4D SM-3MG4E SM-3MG4F SM-3MG4G SM-3MG4H SM-3MG4J 
#        3        1        1        2        1        3        2        2 
# SM-3MG4K SM-3MG4L SM-3MG4M SM-3MG4N SM-3MG4O SM-3MG4P SM-3MG4Q SM-3MG4R 
#        3        2        2        2        1        2        2        1 
# SM-3MG4S SM-3MG4T SM-3MG4U SM-3MG4V SM-3MG4W SM-3MG4X SM-3MG4Y SM-3MG51 
#        1        1        2        2        2        1        1        4 
# SM-3MG52 SM-3MG53 SM-3MG54 SM-3MG55 SM-3MG56 SM-3MG57 SM-3MG58 SM-3MG59 
#        1        1        1        2        1        3        1        2 
# SM-3MG5A SM-3MG5B SM-3MG5C SM-3MG5D SM-3MG5E SM-3MG5F SM-3MG5G SM-3MG5H 
#        1        2        3        3        1        1        2        2 
# SM-3MG5J SM-3MG5L SM-3MG5M SM-3MG5N SM-3MG5O SM-3MG5P SM-3MG5Q SM-3MG5R 
#        2        2        2        2        1        1        1        1 
# SM-3MG5S SM-3MG5T SM-3MG5U SM-3MG5V SM-3MG5W SM-3MG5Y SM-3MG5Z SM-3MG61 
#        1        1        1        1        1        1        1        1 
# SM-3MG62 SM-3MG63 SM-3MG64 SM-3MG65 SM-3MG66 SM-3MG67 SM-3MG68 SM-3MG69 
#        1        1        2        2        1        1        2        1 
# SM-3MG6A SM-3MG6B SM-3MGPL SM-3MGPN SM-3MGPO SM-3MGPP SM-3MGPQ SM-3MGPR 
#        1        2        2        1        2        1        1        1 
# SM-3MGPS SM-3MGPU SM-3MGPV SM-3MGPW SM-3MGPX SM-3MGPY SM-3MGPZ 
#        1        1        1        1        1        1        1 

#How many reads file each subjects has?
with(design,table(tapply(FileName, INDEX=Subject, 
    function(x){length(unique(x))})))
#  8 10 12 14 16 
#  7 20 34 25  9
with(design, tapply(FileName, INDEX=Subject, function(x){length(unique(x))}))
# SM-3MG3L SM-3MG3M SM-3MG3N SM-3MG3O SM-3MG3P SM-3MG3R SM-3MG3U SM-3MG3V 
#       14        8        8       10       12       10       14       14 
# SM-3MG3Y SM-3MG3Z SM-3MG45 SM-3MG46 SM-3MG47 SM-3MG48 SM-3MG49 SM-3MG4A 
#       10       16       16       14       12       14        8       10 
# SM-3MG4B SM-3MG4C SM-3MG4D SM-3MG4E SM-3MG4F SM-3MG4G SM-3MG4H SM-3MG4J 
#       16        8       12       14       12       12       10       10 
# SM-3MG4K SM-3MG4L SM-3MG4M SM-3MG4N SM-3MG4O SM-3MG4P SM-3MG4Q SM-3MG4R 
#       14       10       12       14       12       14       14       12 
# SM-3MG4S SM-3MG4T SM-3MG4U SM-3MG4V SM-3MG4W SM-3MG4X SM-3MG4Y SM-3MG51 
#       10       10       12       14       12       14       12       14 
# SM-3MG52 SM-3MG53 SM-3MG54 SM-3MG55 SM-3MG56 SM-3MG57 SM-3MG58 SM-3MG59 
#        8       14       14       12       10       16       14       12 
# SM-3MG5A SM-3MG5B SM-3MG5C SM-3MG5D SM-3MG5E SM-3MG5F SM-3MG5G SM-3MG5H 
#       10       14       16       16       10       12       14       12 
# SM-3MG5J SM-3MG5L SM-3MG5M SM-3MG5N SM-3MG5O SM-3MG5P SM-3MG5Q SM-3MG5R 
#       10       12       12       14       12       12       12       10 
# SM-3MG5S SM-3MG5T SM-3MG5U SM-3MG5V SM-3MG5W SM-3MG5Y SM-3MG5Z SM-3MG61 
#       10        8       14       14       14       14       12       12 
# SM-3MG62 SM-3MG63 SM-3MG64 SM-3MG65 SM-3MG66 SM-3MG67 SM-3MG68 SM-3MG69 
#       12       10       12       12       12       12       16       10 
# SM-3MG6A SM-3MG6B SM-3MGPL SM-3MGPN SM-3MGPO SM-3MGPP SM-3MGPQ SM-3MGPR 
#       14       16       16        8       10       12       14       12 
# SM-3MGPS SM-3MGPU SM-3MGPV SM-3MGPW SM-3MGPX SM-3MGPY SM-3MGPZ 
#       12       10       12       12       12       12       10
 
################################################################################
##Adding the sequencing metadata form the fq.gz themselves
################################################################################
##How many read sequence per File??
##How many total reads each subject has sequence per File??
##All the read sequenced are OK in the raw data file??

