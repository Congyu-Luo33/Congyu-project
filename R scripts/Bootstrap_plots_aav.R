generate.bootstrap.plots.for.transcriptome.a <- function(sct_data,tt,file_name,time_point,thresh=250,annotLevel=1,reps,full_results=NA,listFileName="",showGNameThresh=25){
  args <- commandArgs(trailingOnly = TRUE)
  # Check the arguments
  correct_length = length(full_results)==5
  required_names = c("joint_results","hit.cells.up","hit.cells.down","bootstrap_data.up","bootstrap_data.down")
  all_required_names = sum(names(full_results) %in% required_names)==5
  if(!correct_length | !all_required_names){stop("ERROR: full_results is not valid output from the ewce_expression_data function. This function only takes data generated from transcriptome analyses.")}
  if(thresh>(2*dim(tt)[1])){stop("ERROR: threshold length must not be more than twice the length of the top table")}
  # Check top table has an MGI.symbol column
  if(!sum(colnames(tt)=="MGI.symbol")==1){stop("ERROR: top table (tt) must have 'MGI.symbol' column")}
  
  for(dirS in c("Up","Down")){
    a = full_results$joint_results
    results = a[as.character(a$Direction)==dirS,]
    
    # Drop genes lacking expression data
    if(dirS=="Up"){tt=tt[order(tt$t,decreasing=TRUE),]}
    if(dirS=="Down"){tt=tt[order(tt$t,decreasing=FALSE),]}
    mouse.hits = as.character(unique(tt$MGI.symbol[1:thresh]))
    mouse.hits = mouse.hits[mouse.hits %in% rownames(sct_data[[1]]$specificity)]
    mouse.bg = as.character(unique(tt$MGI.symbol))
    mouse.bg = mouse.bg[!mouse.bg %in% mouse.hits]
    mouse.bg = mouse.bg[mouse.bg %in% rownames(sct_data[[1]]$specificity)]
    combinedGenes = unique(c(mouse.hits, mouse.bg))
    
    # Get expression data of bootstrapped genes
    signif_res = as.character(results$CellType)
    nReps = 1000
    exp_mats = list()
    #for(cc in signif_res){
    for(cc in as.character(unique(a$CellType))){
      
      exp_mats[[cc]] = matrix(0,nrow=nReps,ncol=length(mouse.hits))
      rownames(exp_mats[[cc]]) = sprintf("Rep%s",1:nReps)
    }
    
    for(s in 1:nReps){
      bootstrap_set = sample(combinedGenes,length(mouse.hits))
      ValidGenes = rownames(sct_data[[annotLevel]]$specificity)[rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set]
      
      expD = sct_data[[annotLevel]]$specificity[ValidGenes,]
      
      #for(cc in signif_res){
      for(cc in as.character(unique(a$CellType))){
        exp_mats[[cc]][s,] = sort(expD[,cc])
      }
    }
    
    # Get expression levels of the hit genes
    hit.exp = sct_data[[annotLevel]]$specificity[mouse.hits,]	#cell.list.exp(mouse.hits)
    
    #print(hit.exp)
    
    graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
      theme(panel.grid.major = element_line(size = .5, color = "grey"),
            axis.line = element_line(size=.7, color = "black"),legend.position = c(0.75, 0.7), text = element_text(size=18),
            axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6)) + theme(legend.title=element_blank())#+ylim(c(1,100))
    
    if (!file.exists("BootstrapPlots")){
      dir.create(file.path(getwd(), "BootstrapPlots"))
    }
    
    tag = sprintf("thresh%s__dir%s",thresh,dirS)
    
    # Plot the QQ plots
    for(cc in as.character(unique(a$CellType))){
      #cc = "Vascular and Leptomeningeal Cells"
      #mean_boot_exp = apply(exp_mats[[cc]],2,mean)
      #mean_boot_exp
      #cc = "Pericytes"
      mean_boot_exp = apply(exp_mats[[cc]],2,mean)
      #mean_boot_exp
      
      
      hit_exp = sort(hit.exp[,cc])
      hit_exp_names = rownames(hit.exp)[order(hit.exp[,cc])]#names(hit_exp)
      dat = data.frame(boot=mean_boot_exp,hit=hit_exp,Gnames=hit_exp_names)
      dat$hit = dat$hit*100
      dat$boot = dat$boot*100
      maxHit = max(dat$hit)
      maxX = max(dat$boot)+0.1*max(dat$boot)
      
      basic_graph = ggplot(dat,aes_string(x="boot",y="hit"))+geom_point(size=1)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")
      
      
      # Plot without text
      pdf(sprintf("BootstrapPlots_aav/%s/%s_%s/qqplot_noText_%s____%s____%s.pdf",time_point,time_point,file_name,tag,listFileName,cc),width=10,height=8)
      print(basic_graph+ggtitle(cc))
      dev.off()
      
      # If a gene has over 25% of it's expression proportion in a celltype, then list the genename
      dat$symLab = ifelse(dat$hit>showGNameThresh,sprintf("  %s", dat$Gnames),'')
      
      #basic_graph = ggplot(dat,aes(x=boot,y=hit))+geom_point(size=2)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
      basic_graph = ggplot(dat,aes_string(x="boot",y="hit"))+geom_point(size=2)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
        geom_abline(intercept = 0, slope = 1, colour = "red")
      
      # Plot with bootstrap distribution
      melt_boot = melt(exp_mats[[cc]])
      colnames(melt_boot) = c("Rep","Pos","Exp")
      actVals = data.frame(pos=as.factor(1:length(hit_exp)),vals=hit_exp)
      
      # Plot with LOG bootstrap distribution
      # - First get the ordered gene names
      rownames(dat)=dat$Gnames
      datOrdered = data.frame(GSym=rownames(dat),Pos=1:dim(dat)[1])
      
      # - Arrange the data frame for plotting
      melt_boot = melt(exp_mats[[cc]])
      colnames(melt_boot) = c("Rep","Pos","Exp")
      melt_boot$Exp = melt_boot$Exp*100
      melt_boot = merge(melt_boot,datOrdered,by="Pos")
      melt_boot$GSym = factor(as.character(melt_boot$GSym),levels=as.character(datOrdered$GSym))
      
      # - Prepare the values of the list genes to be plotted as red dots
      actVals = data.frame(Pos=as.factor(1:length(hit_exp)),vals=hit_exp*100)
      actVals = merge(actVals,datOrdered,by="Pos")
      actVals$GSym = factor(as.character(actVals$GSym),levels=as.character(datOrdered$GSym))
      
      # - Determine whether changes are significant
      p = rep(1,max(melt_boot$Pos))
      for(i in 1:max(melt_boot$Pos)){
        p[i] = sum(actVals[actVals$Pos==i,"vals"]<melt_boot[melt_boot$Pos==i,"Exp"])/length(melt_boot[melt_boot$Pos==i,"Exp"])
      }
      ast = rep("*",max(melt_boot$Pos))
      ast[p>0.05] = ""
      actVals = cbind(actVals[order(actVals$Pos),],ast)
      # - Plot the graph!
      wd = 1+length(unique(melt_boot[,4]))*0.2
      pdf(sprintf("BootstrapPlots_aav/%s/%s_%s/bootDists_LOG_%s___%s____%s.pdf",time_point,time_point,file_name,tag,listFileName,cc),width=wd,height=4)
      #png(sprintf("BootstrapPlots/bootDists_LOG_%s___%s____%s.png",tag,listFileName,cc),width=wd*100,height=4*100)
      melt_boot$GSym = factor(melt_boot$GSym,levels=rev(levels(melt_boot$GSym)))
      print(ggplot(melt_boot)+geom_boxplot(aes_string(x="GSym",y="Exp"),outlier.size=0)+graph_theme+
              theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x=element_text(hjust=0))+
              geom_point(aes_string(x="GSym",y="vals"),col="red",data=actVals)+
              geom_text(aes_string(x="GSym",y="vals",label="ast"),colour="black",col="black",data=actVals)+
              ylab("Cell type specificity\n")+
              xlab("Most specific --> Least specific")+scale_y_log10(limits=c(1,100)))
      
      dev.off()
    }
  }
}
