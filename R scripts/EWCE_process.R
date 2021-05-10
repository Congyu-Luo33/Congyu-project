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
#EWCE
load(file='ctd.rda')
EWCE_g93a <- list.files(path='Results_bootstrap/Tables/p120/', pattern='*.csv')

#This for loop is just to read different files
for (i in EWCE_g93a){
  #i= p70_Vent_Med_White.csv
  t_data_EWCE_test <- read.csv(paste('Results_bootstrap/Tables/p120/', i, sep=''))
  EWCE_results_test <- ewce_expression_data(sct_data=ctd,tt=t_data_EWCE_test, thresh=250, 
                                            reps=10000, annotLevel=1, ttSpecies="mouse",sctSpecies="mouse")
  EWCE_results_test$joint_results$p.adj = p.adjust(EWCE_results_test$joint_results$p)
  save(EWCE_results_test, file=sprintf('test', paste('EWCE_', gsub('.csv','', gsub('DE.','',i)), '.rda', sep='')))
  EWCE_results_test$joint_results[order(EWCE_results_test$joint_results$sd_from_mean),]
  write.csv(EWCE_results_test$joint_results, file=sprintf('Results_wt_g93a/Tables/p120/%s', paste(gsub('DE.','',i), sep='')))
  pdf(file=paste("Results_wt_g93a/Figures/p120/", 'EWCE_', gsub('.csv','',i), '.pdf', sep=''))
  print(ewce.plot(EWCE_results_test$joint_results))
  
  dev.off()
}
