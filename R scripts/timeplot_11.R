library(R.utils)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra)
library(plyr)

#collaborated with get_file_name.sh
#setwd('/Users/apple/Desktop/EWCE')
mkdirs('/Users/apple/Desktop/EWCE/Results_wt_g93a/Figures/timeplot_11')

color10 <- c("#6DC7BE", "#7F3F97", "#2E8B57", "#808184", "#BBBDBF", "#A0522D", "#524FA0", "#EC1C24", "#FFD700", "#FAAF40")
vascular <- c('Pericytes', 'Vascular and Leptomeningeal Cells', 'Vascular Endothelial', 'Vascular Smooth Muscle Cell')
parenchymal <- c('Astrocytes', 'Interneurons', 'Microglia', 'Oligodendrocyte Precursor', 'Oligodendrocytes', 'Pyramidal Neurons')


args <- commandArgs(trailingOnly = TRUE)


fileDir <- '/Users/apple/Desktop/EWCE/Results_wt_g93a/Tables/'
aav30 <- read.csv(sprintf('%sp30/p30_%s.csv',fileDir,args[1]),header=TRUE)
aav30$timepoint <- as.numeric('30')
aav70 <- read.csv(sprintf('%sp70/p70_%s.csv',fileDir,args[1]),header=TRUE)
aav70$timepoint <- as.numeric('70')
aav100 <- read.csv(sprintf('%sp100/p100_%s.csv',fileDir,args[1]),header=TRUE)
aav100$timepoint <- as.numeric('100')
aav120 <- read.csv(sprintf('%sp120/p120_%s.csv',fileDir,args[1]),header=TRUE)
aav120$timepoint <- as.numeric('120')

data_aav <- rbind(aav30, aav70, aav100, aav120)
data_aav$sd_from_mean[data_aav$sd_from_mean < 0] <- 0
data_aav$significance[data_aav$p.adj <= 0.05] <- '*'
data_aav$significance[data_aav$p.adj > 0.05] <- ''
data_aav$category[data_aav$CellType %in% vascular] <- 'vascular'
data_aav$category[data_aav$CellType %in% parenchymal] <- 'parenchymal'
data_aav$category <- relevel(factor(data_aav$category), ref='vascular')
data_aav$sd_from_mean[data_aav$Direction=='Down'] <- data_aav$sd_from_mean[data_aav$Direction=='Down']*(-1)
data_aav$Direction <- factor(data_aav$Direction, levels=c('Up', 'Down'))

direction_wt2 <- data_aav$Direction
for(j in direction_wt2){
  plot_up <- data_aav[data_aav$Direction=='Up',]
  plot_down <- data_aav[data_aav$Direction=='Down',]
  if(j=='Up'){
    p1=ggplot(data=plot_up, aes(x=factor(timepoint), y=sd_from_mean, color=CellType, group=CellType)) + 
      geom_point() +  geom_line() + scale_color_manual(values = color10) + 
      geom_hline(yintercept = 0, color='grey') +
      geom_text(label=plot_up$significance, nudge_y = 0.5, size=10, color='black')+
      facet_grid(scales='free_y', space = 'free_y') + xlab(NULL)  + ylab(NULL) +
      scale_y_continuous(limits = c(0,40), breaks=c(0,10,20,30,40))+theme_bw()+
      ggtitle(sprintf('%s',args[1])) +
      theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5, face = 'bold', size = rel(1.3)), strip.background.y = element_blank(),
            strip.text.y = element_blank(),panel.border=element_rect(color="white"),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=16,face = "bold"))
  }
  else{
    p2=ggplot(data=plot_down, aes(x=factor(timepoint), y=sd_from_mean, color=CellType, group=CellType)) +
      geom_point() + geom_line() + scale_color_manual(values = color10) + 
      geom_hline(yintercept = 0, color='grey') +
      geom_text(label=plot_down$significance, nudge_y = 0.5, size=10, color='black')+
      facet_grid(scales='free_y', space = 'free_y')  + ylab(NULL) + xlab(NULL) +
      scale_y_continuous(limits=c(-40,0),breaks=c(-40,-30,-20,-10,0))+theme_bw()+
      theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5, face = 'bold', size = rel(1.3)), strip.background.y = element_blank(),
            strip.text.y = element_blank(),panel.border=element_rect(color="white"),axis.text=element_text(size=16,face = "bold"))
    
  }
}
pdf(file=sprintf('/Users/apple/Desktop/EWCE/Results_wt_g93a/Figures/timeplot_11_2/%s_timeplot.pdf',args[1]), width=5)

plot <- plot_grid(p1, p2, ncol = 1,nrow=2,align = 'v')

y.grob <- textGrob('z-score', 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

x.grob <- textGrob('Days',
                   gp=gpar(fontface="bold", col="black", fontsize=15),hjust = 2.5)

print(grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob)))

dev.off()