
# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : notintable
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : check which genes from TregvsTconv_Thy_DEGnoco.txt
# do not appear in the heatmap
#
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      first prototype
#
#
# **********************************************************************




notintable<-function(AllClusters){
  
  #deg datasets
  genesTreg <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  
  genes<- genesTreg$hgnc_symbol
  genes<-toupper(genes)
  genes<-unique(genes)
  
  GenesSet<-AllClusters$hgnc_symbol
  GenesSet<-toupper(GenesSet)
  GenesSet<-unique(GenesSet)
  
  Genes0s<-setdiff(genes,GenesSet)
  return(Genes0s)
  
}