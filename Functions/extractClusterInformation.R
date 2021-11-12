

# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : Extract cluster info
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Extract cluster info from the heatmap results
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      first prototype
#
#
# **********************************************************************




#diffbinding
extractClusterInformation_DiffBind<-function(HM){
  #___________ cluster Genes
  
  r.dend <- row_dend(HM)  #If needed, extract row dendrogram
  rcl.list <- row_order(HM)  #Extract clusters (output is a list)
  
  lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters
  
  
  library(magrittr) # needed to load the pipe function '%>%'
  
  clu_df <- lapply(names(rcl.list), function(i){
    r=rownames(data_wideHeatmap)
    out <- data.frame(ensembl_gene_id = r[rcl.list[[i]]],
                      ClusterGene = paste0(i),
                      stringsAsFactors = FALSE)
    return(out)
  }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
    do.call(rbind, .)
  
  ClustersGenes<-unique(data.frame(clu_df))
  
  #_______
  
  #___________ cluster TFBS
  
  c.dend <- column_dend(HM)  #If needed, extract row dendrogram
  cl.list <- column_order(HM)  #Extract clusters (output is a list)
  
  lapply(cl.list, function(x) length(x))  #check/confirm size gene clusters
  
  library(magrittr) # needed to load the pipe function '%>%'
  
  clu_TFBS<- lapply(names(cl.list), function(i){
    t=colnames(data_wideHeatmap)
    out <- data.frame(TFBS_name = t[cl.list[[i]]],
                      ClusterTFBS = paste0(i),
                      stringsAsFactors = FALSE)
    return(out)
  }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
    do.call(rbind, .)
  
  
  ClustersTFBS<-unique(as.data.frame(clu_TFBS))
  
  
  
  All_clusters<-merge(data, ClustersGenes, by="ensembl_gene_id")
  All_clusters<-as.data.frame(All_clusters)
  
  #All_clusters$TFBS_name <- gsub(" _ ", "_", All_clusters$TFBS_name)
  
  
  
  All_clusters<-inner_join(All_clusters,ClustersTFBS,by="TFBS_name")
  
  
  #obtain hgnc_symbol
  obtainSymbol <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  
  
  obtainSymbol<-data.frame(obtainSymbol$ensembl_gene_id, obtainSymbol$hgnc_symbol)
  names(obtainSymbol)<-c("ensembl_gene_id", "hgnc_symbol")
  
  
  All_clusters<-merge(x = All_clusters, y = obtainSymbol, by = "ensembl_gene_id", all.x = TRUE)
  
  
  
  All_clusters<-data.frame(All_clusters$ensembl_gene_id, All_clusters$hgnc_symbol, All_clusters$TFBS_name, All_clusters$meanDiffBind, All_clusters$ClusterGene, All_clusters$ClusterTFBS)
  
  names(All_clusters)<-c("ensembl_gene_id", "hgnc_symbol", "TFBS_name", "meanDiffBind", "ClusterGene", "ClusterTFBS")
  
  
  
  All_clusters<-unique(All_clusters)
  
  
  
  #divideClusters
  
  CLustersGenesData<-data.frame(All_clusters$ensembl_gene_id, All_clusters$hgnc_symbol, All_clusters$ClusterGene)
  
  names(CLustersGenesData)<-c("ensembl_gene_id", "hgnc_symbol", "ClusterGene")
  
  CLustersGenesData<<-unique(CLustersGenesData)
  
  ClustersTFBS<-data.frame(All_clusters$TFBS_name, All_clusters$ClusterTFBS)
  
  names(ClustersTFBS)<-c("TFBS_name", "ClusterTFBS")
  
  ClustersTFBS<<-unique(ClustersTFBS)
  
  
  
  returnlist<-list(All_clusters=All_clusters, CLustersGenesData=CLustersGenesData, ClustersTFBS=ClustersTFBS)
  
  return(returnlist)
  
}




#diffbinding - hgnc
extractClusterInformation_DiffBind_HGNC<-function(HM){
  #___________ cluster Genes
  
  r.dend <- row_dend(HM)  #If needed, extract row dendrogram
  rcl.list <- row_order(HM)  #Extract clusters (output is a list)
  
  lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters
  
  library(magrittr) # needed to load the pipe function '%>%'
  
  clu_df <- lapply(names(rcl.list), function(i){
    r=rownames(data_wideHeatmap)
    out <- data.frame(hgnc_symbol = r[rcl.list[[i]]],
                      ClusterGene = paste0(i),
                      stringsAsFactors = FALSE)
    return(out)
  }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
    do.call(rbind, .)
  
  ClustersGenes<-unique(data.frame(clu_df))
  
  #_______
  
  #___________ cluster TFBS
  
  c.dend <- column_dend(HM)  #If needed, extract row dendrogram
  cl.list <- column_order(HM)  #Extract clusters (output is a list)
  
  lapply(cl.list, function(x) length(x))  #check/confirm size gene clusters
  
  library(magrittr) # needed to load the pipe function '%>%'
  
  clu_TFBS<- lapply(names(cl.list), function(i){
    t=colnames(data_wideHeatmap)
    out <- data.frame(TFBS_name = t[cl.list[[i]]],
                      ClusterTFBS = paste0(i),
                      stringsAsFactors = FALSE)
    return(out)
  }) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
    do.call(rbind, .)
  
  
  ClustersTFBS<-unique(as.data.frame(clu_TFBS))
  
  
  All_clusters<-merge(data, ClustersGenes, by="hgnc_symbol")
  All_clusters<-as.data.frame(All_clusters)
  
  #All_clusters$TFBS_name <- gsub(" _ ", "_", All_clusters$TFBS_name)
  
  
  
  All_clusters<-inner_join(All_clusters,ClustersTFBS,by="TFBS_name")
  
  
  # #obtain hgnc_symbol
  # obtainSymbol <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  # 
  # 
  # obtainSymbol<-data.frame(obtainSymbol$ensembl_gene_id, obtainSymbol$hgnc_symbol)
  # names(obtainSymbol)<-c("ensembl_gene_id", "hgnc_symbol")
  # 
  # 
  # All_clusters<-merge(x = All_clusters, y = obtainSymbol, by = "ensembl_gene_id", all.x = TRUE)
  # 
  
  
  # All_clusters<-data.frame(All_clusters$ensembl_gene_id, All_clusters$hgnc_symbol, All_clusters$TFBS_name, All_clusters$meanDiffBind, All_clusters$ClusterGene, All_clusters$ClusterTFBS)
  # 
  # names(All_clusters)<-c("ensembl_gene_id", "hgnc_symbol", "TFBS_name", "meanDiffBind", "ClusterGene", "ClusterTFBS")
  # 
  # 
  # 
  All_clusters<-unique(All_clusters)
  
  
  
  #divideClusters
  
  CLustersGenesData<-data.frame(All_clusters$ensembl_gene_id, All_clusters$hgnc_symbol, All_clusters$ClusterGene)
  
  names(CLustersGenesData)<-c("ensembl_gene_id", "hgnc_symbol", "ClusterGene")
  
  CLustersGenesData<<-unique(CLustersGenesData)
  
  ClustersTFBS<-data.frame(All_clusters$TFBS_name, All_clusters$ClusterTFBS)
  
  names(ClustersTFBS)<-c("TFBS_name", "ClusterTFBS")
  
  ClustersTFBS<<-unique(ClustersTFBS)
  
  
  
  returnlist<-list(All_clusters=All_clusters, CLustersGenesData=CLustersGenesData, ClustersTFBS=ClustersTFBS)
  
  return(returnlist)
  
}