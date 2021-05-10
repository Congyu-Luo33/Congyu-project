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
library(MASS) 
library(reshape2) 
library(reshape) 


args <- commandArgs(trailingOnly = TRUE)
time <- as.character('p120')
mkdirs(sprintf("BootstrapPlots_aav/%s/%s_%s",time,time,args[1]))
#collaborated with get_file_name_for_bootstrap.sh
#bootstrap_plot
load(file='ctd.rda')
source('ewce_expression_data_aav.r')
source('Bootstrap_plots_aav.R')

t1<-read.csv(sprintf('Results_bootstrap/Tables/%s/%s_%s.csv',time,time,args[1]),header=TRUE)
t2<-ewce_expression_data(sct_data=ctd,tt=t1)
generate.bootstrap.plots.for.transcriptome.a(sct_data=ctd,tt=t1,full_results = t2,
                                             time_point = time,file_name=args[1])
