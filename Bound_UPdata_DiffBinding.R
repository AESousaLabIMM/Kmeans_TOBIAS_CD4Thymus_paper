# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : Bound_UPdata_DiffBinding 
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Code to Create the kmeans heatmap for diffbinding in 
# the UP DEG data 
#
#
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110   Susana      1      first prototype
#
#
# **********************************************************************



# **********************************************************************
# Packages and Functions
# **********************************************************************

library (plyr)
library(readr)
library(ggplot2)
library(tibble)
library(tidyverse)
library(circlize)
library(tidyr)
library(magrittr)

#source("Functions/CreateMatrix.R")
source("Functions/extractClusterInformation.R")
source("Functions/notintable.R")
source("Functions/ExtractInfoFromDataset.R")
source("Functions/AssessNumClusters.R")
source("Functions/ColoursHeatmap.R")

# **********************************************************************
# Extract and Process Data
# **********************************************************************

#### -- The input is the dataset extracted from the TOBIAS framework for up regulated genes
# just replace it in the code below


#extract data and subset for treg_bound==1
data <- read.table("TobiasOutputTable_UpregulatedGenes.txt", sep = "\t", header = TRUE )
data <- subset(data, treg_bound==1)


#run function to process data
returnlist<-MatrixAndBar_DiffBinding(data)

#output data
data<-returnlist[[1]]

#output matrix for heatmap
data_wide<-returnlist[[2]]

#output expression data for rows - logFC
RowBar<-returnlist[[3]]

#output expression data for cols - logFC
ColData<-returnlist[[4]]




# **********************************************************************
# Calculate ideal number of clusters for Row and Col
# **********************************************************************

AssessNumClusters(data_wide)


# **********************************************************************
# set ideal number of clusters 
# **********************************************************************
GeneClusters<-6
TFBSClusters<-7


# **********************************************************************
#  Plot Heatmap with clustering - scaling data by rows
# **********************************************************************

#necessary packages
library(cluster)
library(ComplexHeatmap)
library(scales)
library(colorRamps)
library(RColorBrewer)
options(scipen = 7)





#scale data by row
data_wideHeatmap<-t(scale(t(data_wide)))


# get colour scheme
mycols<-getcolours_diffBindingHeatmap_RowScaling(data_wideHeatmap)



#row annotations - bar plot
ha2 = rowAnnotation(
  ExpGene = anno_barplot(
    RowBar$expGene, 
    bar_width = 1, 
    width = unit(3, "cm"),
    gp = gpar(fill = RowBar$color, col = RowBar$color),
    border = FALSE
  ), show_annotation_name = TRUE)


#col annotation - bar plot
ha1 = HeatmapAnnotation(
  ExpTF = anno_barplot(
    ColData$expTF, 
    bar_width = 1, 
    height = unit(2, "cm")
  ), show_annotation_name = TRUE)





# Annotations


#row annotations - bar plot
ha2 = rowAnnotation(
  ExpGene = anno_barplot(
    RowBar$expGene, 
    bar_width = 1, 
    width = unit(3, "cm"),
    gp = gpar(fill = RowBar$color, col = RowBar$color),
    border = FALSE
  ), show_annotation_name = TRUE)


#col annotation - bar plot
ha1 = HeatmapAnnotation(
  ExpTF = anno_barplot(
    ColData$expTF, 
    bar_width = 1, 
    height = unit(2, "cm")
  ), show_annotation_name = TRUE)



# ___________ create heatmap __________

set.seed(123)

HMr<-ha2+Heatmap(data_wideHeatmap,
                 name = "Bound_UPdata_DiffBinding_RowScaling",
                 col = mycols,
                 row_dend_reorder=TRUE,
                 row_km = GeneClusters,
                 column_km = TFBSClusters,
                 row_km_repeats = 100,
                 column_km_repeats = 100,
                 top_annotation = ha1,
                 show_row_names = FALSE
)
#__________________________________________________

#___plot in svg

svg("Bound_UPdata_DiffBinding_RowScaling.svg",h=20, w=30)

HMr<-draw(HMr)

dev.off()


#___plot in pdf

pdf("Bound_UPdata_DiffBinding_RowScaling.pdf",h=35, w=35)


HMr<-draw(HMr)

dev.off()


#___plot in png

png("Bound_UPdata_DiffBinding_RowScaling.png",h=3500, w=3500)

HMr<-draw(HMr)

dev.off()



# **********************************************************************
#  Extract clustering data from the heatmap
# **********************************************************************

#extract clustering data 
returnlistclusters<-extractClusterInformation_DiffBind(HMr)

#____all clusters data
AllClusters<-returnlistclusters[[1]]

write.csv(AllClusters, file="Bound_UPdata_DiffBinding_RowScaling_AllClusters.csv")

#____gene clusters data

ClusterGenes<-returnlistclusters[[2]]

write.csv(ClusterGenes, file="Bound_UPdata_DiffBinding_RowScaling_ClustersGenes.csv")



#____TFBS clusters data

ClustersTFBS<-returnlistclusters[[3]]

write.csv(ClustersTFBS, file="Bound_UPdata_DiffBinding_RowScaling_ClustersTFSB.csv")






# **********************************************************************
#  Check which genes from the List cannot be found in the heatmap 
# **********************************************************************
Genes0s<-notintable(AllClusters)

write.csv(Genes0s, "Bound_UPdata_DiffBinding_RowScaling_Genes0s.csv")

