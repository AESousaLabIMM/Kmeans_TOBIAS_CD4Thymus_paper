# **********************************************************************
# Project           : Kmeans_TOBIAS_CD4Thymus
#
# Program name      : ColoursHeatmap
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Get the colours for the heatmap
#
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      first prototype
#
#
# **********************************************************************



# ----- Function to set colours for Diff Binding



getcolours_diffBindingHeatmap_RowScaling<-function(datamatrix){
  
  Modes<-VanillaICE::rowModes(datamatrix)
  Minimum<-min(Modes)
  Maximum<-max(Modes)
  
  
  library(circlize)
  col_fun = colorRamp2(c(min(datamatrix)-.001,Minimum-0.1, Minimum-0.01, Minimum:Maximum, Maximum+0.01, Maximum+0.1,max(datamatrix)+.001),
                       c("#4fff2e", "#189b00", "#106a00", "black", "black" ,"#b36200", "#e67e00", "#ffa333"), space = "XYZ") #green black red
  
  
  
  return(col_fun)
}



