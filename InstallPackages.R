# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : InstallPackages 
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Necessary Packages to Run this Project
#
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110   Susana      1      first prototype
#
#
# **********************************************************************

#basic packages
install.packages('knitr')
install.packages('plyr')
install.packages('readr')
install.packages('ggplot2')
install.packages('tibble')
install.packages('tidyverse')
install.packages('colortools')
install.packages('tidyr')
install.packages('magrittr')
install.packages('gtools')

install.packages('factoextra')
install.packages('circlize')


install.packages('cluster')
install.packages('scales')
install.packages('RColorBrewer')


install.packages('colortools')
install.packages('colorRamps')



#ComplexHeatmap
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


# BSgenome.Hsapiens.UCSC.hg18

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg18")


# VanillaICE
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VanillaICE")




#oligoClasses
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligoClasses")



#fgsea

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
