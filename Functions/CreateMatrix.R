
# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : CreateMAtrix
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Extract data from the original tobias dataset and create
# - table with data 
# - table for heatmap
# - row expression data
# - column expression data
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      first prototype
#
#
# **********************************************************************

############### DiffBind


# extract matrix for heatmap from tobias data - data in matrix is DiffBind
CreateMatrix_DiffBind<-function(data){
  
  #necessary packages
  library (plyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(tidyverse)
  library(tidyr)
  library(magrittr)
  
  #All of the Genes 
  genesTreg <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  
  #_____processing
  
  #correct ensemble
  data$ensembl_gene_id<-data$gene_id
  
  #remove rows where score is 0
  data=subset(data,data$tconv_score!=0)
  
  #normalize names
  if(!"TFBS_motif" %in% colnames(data))
  {
    x<-c("TFBS_motif")
    data[x[!(x %in% colnames(data))]] = data$TFBS_name
  }
  
  #just genes 
  DEGS<-data.frame(genesTreg$ensembl_gene_id)
  names(DEGS)<-c("ensembl_gene_id")
  data<-merge(data, DEGS, by= "ensembl_gene_id")
  
  #prepare data
  data<-data.frame(data$TFBS_name, data$treg_tconv_log2fc, data$ensembl_gene_id)
  names(data)<-c( "TFBS_name", "DiffBind", "ensembl_gene_id")
  
  #to retrieve duplicates
  data$symbolTFBS <- paste(data$ensembl_gene_id,"_",data$TFBS_name)
  
  
  #reformat data
  data<-data.frame(data$ensembl_gene_id, data$TFBS_name, data$DiffBind, data$symbolTFBS)
  names(data)<-c("ensembl_gene_id", "TFBS_name", "DiffBind", "symbolTFBS")
  
  
  #means
  means<-data.frame(tapply(data$DiffBind, data$symbolTFBS, mean))
  library(dplyr)
  means <- tibble::rownames_to_column(means, "symbolTFBS")
  colnames(means)<-c("symbolTFBS","meanDiffBind")
  
  data<-merge(data, means, by='symbolTFBS')
  
  
  data = data[order(data[,'symbolTFBS'],-data[,'meanDiffBind']),]
  data = data[!duplicated(data$symbolTFBS),]
  
  
  #convert for heatmap
  data<-data.frame(data$ensembl_gene_id, data$TFBS_name, data$meanDiffBind)
  names(data)<-c("ensembl_gene_id", "TFBS_name", "meanDiffBind")
  
  
  #just necessary
  data<-data.frame(data$ensembl_gene_id, data$TFBS_name, data$meanDiffBind)
  names(data)<-c("ensembl_gene_id", "TFBS_name", "meanDiffBind" )
  
  
  data$symbolTFBS <- paste(data$ensembl_gene_id,"_",data$TFBS_name)
  
  data = data[order(data[,'TFBS_name'],-data[,'meanDiffBind']),]
  data = data[!duplicated(data$symbolTFBS),]
  
  #just necessary
  data<-data.frame(data$ensembl_gene_id, data$TFBS_name, data$meanDiffBind)
  names(data)<-c("ensembl_gene_id", "TFBS_name", "meanDiffBind" )
  
  #create matrix
  data_wide <- spread(data, TFBS_name, meanDiffBind)
  data_wide<-as.data.frame(data_wide)
  
  
  data_wide<-data_wide[!is.na(data_wide$ensembl_gene_id),]
  rownames(data_wide)<-data_wide$ensembl_gene_id
  data_wide<-data_wide[,-1] # delete column 1
  
  data_wide[is.na(data_wide)] <- 0
  data_wide<-data.matrix(data_wide)
  
  # ------- gene col and row data
  
  library(gtools)
  
  ensembl_gene_idS<-data.frame(data$ensembl_gene_id)
  names(ensembl_gene_idS)<-c("ensembl_gene_id")
  
  
  #geneids
  #gene_ids <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  gene_ids <- read.table("Data/ExpressionData/treg_and_tconv_filtExpLog.txt", sep = "\t", header = TRUE )
  
  gene_ids<-data.frame(gene_ids$ensembl_gene_id)
  names(gene_ids)<-c("ensembl_gene_id")
  
  #expgene
  #ExpLog <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  ExpLog <- read.table("Data/ExpressionData/treg_and_tconv_filtExpLog.txt", sep = "\t", header = TRUE )
  
  
  #explog
  explog<-data.frame(ExpLog$ensembl_gene_id, ExpLog$tTreg)
  names(explog)<-c("ensembl_gene_id", "expGene")
  RowBar<-merge(ensembl_gene_idS, explog, by= "ensembl_gene_id")
  
  RowBar<-data.frame(RowBar$ensembl_gene_id, RowBar$expGene)
  names(RowBar)<-c("ensembl_gene_id", "expGene")
  
  
  RowBar$expGeneLinear <- logratio2foldchange(RowBar$expGene, base = 2)
  
  RowBar$expGeneLinear   <- ifelse(RowBar$expGene<0 , 2^(RowBar$expGene),logratio2foldchange(RowBar$expGene, base = 2))
  
  
  RowBar<-unique(RowBar)
  
  
  RowBar$static<-2
  
  
  #logFC <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  logFC <- read.table("Data/ExpressionData/treg_and_tconv_filtExpLog.txt", sep = "\t", header = TRUE )
  
  
  logFC<-data.frame(logFC$ensembl_gene_id, logFC$tTreg)
  names(logFC)<-c("ensembl_gene_id", "tTreg")
  
  RowBar<-merge(logFC, RowBar, by= "ensembl_gene_id")
  
  
  RowBar$type<- ifelse(RowBar$tTreg>0, "UpGene", "DownGene")
  RowBar$color<- ifelse(RowBar$tTreg>0, "red", "blue")
  
  
  #expTF
  
  TFData<-data.frame(data$TFBS_name)
  names(TFData)<-c("TFBS_name")
  
  TFData<-TFData %>% separate(TFBS_name, c("TF", "TFBS"), sep = "_")
  
  
  #ExpTF <- read.table("Data/ExpressionData/TregvsTconv_Thy_DEGnoco.txt", sep = "\t", header = TRUE )
  ExpTF <- read.table("Data/ExpressionData/treg_and_tconv_filtExpLog.txt", sep = "\t", header = TRUE )
  
  
  ExpTF<-data.frame(ExpTF$hgnc_symbol, ExpTF$tTreg)
  names(ExpTF)<-c("TF", "expTF")
  
  ColData<-merge(TFData, ExpTF, by= "TF")
  
  ColData$TFBS_name<-paste(ColData$TF,"_",ColData$TFBS)
  
  
  ColData<-data.frame(ColData$TFBS_name, ColData$expTF)
  names(ColData)<-c("TFBS_name", "expTF")
  
  
  ColData$expTFLinear  <- ifelse(ColData$expTF<0 , 2^(ColData$expTF),logratio2foldchange(ColData$expTF, base = 2))
  ColData<-unique(ColData)
  
  returnlist<-list(data=data,data_wide=data_wide, RowBar=RowBar, ColData=ColData)
  
  return(returnlist)
}

