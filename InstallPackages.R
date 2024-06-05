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


#install colortools (deprecated in R)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/colortools/colortools_0.1.5.tar.gz"
install.packages(packageurl, repos=NULL, type="source")





#ComplexHeatmap
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


# BSgenome.Hsapiens.UCSC.hg18

BiocManager::install("BSgenome.Hsapiens.UCSC.hg18")


# VanillaICE

BiocManager::install("VanillaICE", force=TRUE)




#oligoClasses

BiocManager::install("oligoClasses")



#fgsea

BiocManager::install("fgsea")
