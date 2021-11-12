# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : ExtractInfoFromDataset
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Assess Ideal Num Clusters for Dataset
#
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      first 
#
#
# **********************************************************************


AssessNumClusters<-function(datamatrix){
  library(factoextra)  
  
  # __________________________ Cluster by Rows
  data_widematrixRows<-scale(data_wide, center = TRUE)
  data_widematrixRows<-as.data.frame(data_widematrixRows)
  
  #max number clusters
  kmax<-nrow(data_widematrixRows)-1
  
  if(kmax>10){kmax=10}
  
  
  # Elbow method
  ElbowGenes<-fviz_nbclust(data_widematrixRows, kmeans, iter.max=1000, method = "wss", k.max = kmax) +
    labs(subtitle = "Elbow method - Genes")
  ElbowGenes
  
  
  # Silhouette method
  SilhouetteGenes<-fviz_nbclust(data_widematrixRows,kmeans,iter.max=1000, method = "silhouette", k.max = kmax)+
    labs(subtitle = "Silhouette method - Genes")
  SilhouetteGenes
  
  
  # ______________ Cluster by Columns
  data_widematrixColumns<-scale(t(data_wide))
  data_widematrixColumns[is.nan(data_widematrixColumns)] <- 0
  
  
  kmax<-nrow(data_widematrixColumns)-1
  
  if(kmax>10){kmax=10}
  
  # Elbow method
  ElbowTFBS<-fviz_nbclust(data_widematrixColumns, kmeans, iter.max=1000, method = "wss", k.max = kmax) +
    labs(subtitle = "Elbow method - TFBS")
  ElbowTFBS
  
  
  # Silhouette method
  SilhouetteTFBS<-fviz_nbclust(data_widematrixColumns,kmeans,iter.max=100, method = "silhouette", k.max = kmax)+
    labs(subtitle = "Silhouette method - TFBS")
  
  SilhouetteTFBS
  
  
  return(list(ElbowGenes=ElbowGenes,SilhouetteGenes=SilhouetteGenes, ElbowTFBS= ElbowTFBS, SilhouetteTFBS=SilhouetteTFBS ))
  
}