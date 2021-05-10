library(gdata)
library(tidyfst)
library(tidyr)
library(dplyr)
library(edgeR)
library(dplyr)
library(R.utils)
library(EWCE)
library(ggplot2)
library(ewceData)
library(cowplot)
library(limma)
library(readxl)

mkdirs('Results_bootstrap')
mkdirs('Results_bootstrap/Tables/p120')
mkdirs('Results_bootstrap/Figures/p120')

#load .csv file
data_pt <- read.csv(file = 'aav9776_data/p120_WT_g93a.csv', header = TRUE, row.names=NULL, na.strings = c("inf"))
#extract specific columns by Column number or column name
Vent_Med_White <- data_pt[-1,c(1,2:6)]
Vent_Horn <- data_pt[-1,c(1,7:11)]
Vent_Lat_White <- data_pt[-1,c(1,12:16)]
Med_Grey <- data_pt[-1,c(1,17:21)]
Dors_Horn <- data_pt[-1,c(1,22:26)]
Dors_Edge <- data_pt[-1,c(1,27:31)]
Med_Lat_White <- data_pt[-1,c(1,32:36)]
Vent_Edge <- data_pt[-1,c(1,37:41)]
Dors_Med_White <- data_pt[-1,c(1,42:46)]
Cent_Can <- data_pt[-1,c(1,47:51)]
Lat_Edge <- data_pt[-1,c(1,52:56)]

Anatomy_areas <- list(Vent_Med_White, Vent_Horn, Vent_Lat_White, Med_Grey, Dors_Horn, Dors_Edge, Med_Lat_White, 
                      Vent_Edge, Dors_Med_White,Cent_Can, Lat_Edge ) 
#Anatomy_areas <- list(Vent_Lat_White) 
for(j in Anatomy_areas) {
  file_name = colnames(j[2])
  colnames(j) <- c('gene','BF','wt_mean','g93a_mean','wt_sd','g93a_sd')
  #make sure the type of data is numeric
  wt_mean <- as.numeric(unlist(j$wt_mean))
  g93a_mean <- as.numeric(unlist(j$g93a_mean))
  wt_sd <- as.numeric(unlist(j$wt_sd))
  g93a_sd <- as.numeric(unlist(j$g93a_sd))
  bf <- as.numeric(unlist(j$BF))
  Tt  <- (g93a_mean - wt_mean)/(sqrt((g93a_sd^2 + wt_sd^2)/2))
  j_tf <- transform(j,TF = Tt)
  
  
  #Pick out genes that appear more than once
  duplicate_gene_name_tf_j <- j_tf[duplicated(j_tf$gene),]
  #Pick out rows with duplicate genes
  duplicated_gene_tf_j <- j_tf[j_tf$gene %in% duplicate_gene_name_tf_j$gene,]
  #Each gene keeps only the row with the highest TF value
  duplicated_decreasing_tf_j <- arrange_dt(duplicated_gene_tf_j, -duplicated_gene_tf_j$TF)
  duplicated_decreasing_tf_j <- duplicated_decreasing_tf_j[duplicated(duplicated_decreasing_tf_j$gene),]
  unique_tf_j <- j_tf[!(j_tf$gene %in% duplicated_decreasing_tf_j$gene & j_tf$TF %in% duplicated_decreasing_tf_j$TF),]
  #Renamed the row name to the gene name
  rownames(unique_tf_j) <- unique_tf_j[,1]
  ewce_input <- unique_tf_j[,c('gene','TF','BF')]
  colnames(ewce_input) <- c('MGI.symbol','t','B')
  write.csv(ewce_input,file=paste("Results_bootstrap/Tables/p120/", 'p120_', gsub('.csv','',file_name), '.csv', sep=''), quote=F,row.names = F)
}


