#' Bootstrap celltype enrichment test for transcriptome data
#'
#' \code{ewce_expression_data} takes a differential expression table and determines the probability of cell-type enrichment in the up & down regulated genes
#'
#'
#' @param sct_data List generated using \code{\link{generate.celltype.data}}
#' @param tt Differential expression table. Can be output of limma::topTable function. Minimum requirement is that one column stores a metric of increased/decreased expression (i.e. log fold change, t-statistic for differential expression etc) and another contains either HGNC or MGI symbols.
#' @param sortBy Column name of metric in tt which should be used to sort up- from down- regulated genes. Default="t"
#' @param thresh The number of up- and down- regulated genes to be included in each analysis. Dafault=250
#' @param reps Number of random gene lists to generate (default=100 but should be over 10000 for publication quality results)
#' @param annotLevel an integer indicating which level of the annotation to analyse. Default = 1.
#' @param ttSpecies Either 'mouse' or 'human' depending on which species the differential expression table was generated from
#' @param sctSpecies Either 'mouse' or 'human' depending on which species the single cell data was generated from
#' @return A list containing five data frames:
#' \itemize{
#'   \item \code{results}: dataframe in which each row gives the statistics (p-value, fold change and number of standard deviations from the mean) associated with the enrichment of the stated cell type in the gene list. An additional column *Direction* stores whether it the result is from the up or downregulated set.
#'   \item \code{hit.cells.up}: vector containing the summed proportion of expression in each cell type for the target list
#'   \item \code{hit.cells.down}: vector containing the summed proportion of expression in each cell type for the target list#'
#'   \item \code{bootstrap_data.up}: matrix in which each row represents the summed proportion of expression in each cell type for one of the random lists
#'   \item \code{bootstrap_data.down}: matrix in which each row represents the summed proportion of expression in each cell type for one of the random lists
#' }
#' @examples
#' # Load the single cell data
#' data("ctd")
#'
#' # Set the parameters for the analysis
#' reps=100 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
#' annotLevel=1 # <- Use cell level annotations (i.e. Interneurons)
#'
#' # Load the top table
#' data(tt_alzh)
#'
#' tt_results = ewce_expression_data(sct_data=ctd,tt=tt_alzh,annotLevel=1,
#'     ttSpecies="human",sctSpecies="mouse")
#' @export
# @import reshape2
# @import plyr

ewce_expression_data <- function(sct_data,annotLevel=1,tt,sortBy="t",thresh=250,reps=10000,ttSpecies="mouse",sctSpecies="mouse"){

    # Check the arguments
    if(!sortBy %in% colnames(tt)){stop("ERROR: tt does not contain a column with value passed in sortBy argument")}
    if(dim(tt)[1]<(thresh*2)){stop("ERROR: length of table is less than twice the size of threshold")}

    # Check that the top table has correct columns
    if(ttSpecies=="human" & !"HGNC.symbol" %in% colnames(tt)){stop("ERROR: if ttSpecies==human then there must be an HGNC.symbol column")}
    if(ttSpecies=="mouse" & !"MGI.symbol" %in% colnames(tt)){stop("ERROR: if ttSpecies==mouse then there must be an MGI.symbol column")}

    if(ttSpecies=="human"){
        tt$MGI.symbol = tt$HGNC.symbol
    }
    tt$MGI.symbol = as.character(tt$MGI.symbol)

    tt2=tt

     if(ttSpecies=="human"){
         if(!"HGNC.symbol" %in% colnames(tt)){stop("ERROR: tt does not contain an 'HGNC.symbol' column")}
         # Load the ortholog data
         mouse_to_human_homologs=NULL # <-- THIS LINE ONLY INCLUDED TO FOOL devtools::check()
         data("mouse_to_human_homologs", envir = environment())
         m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])

         # Drop genes without orthologs
         tt2 = merge(tt,m2h,by="HGNC.symbol") # Obtain MGI orthologs
     }else{
         if(!"MGI.symbol" %in% colnames(tt)){stop("ERROR: tt does not contain an 'MGI.symbol' column")}
         tt2 = tt
     }

    # Sort from down-->up regulated
    tt3 = tt2[order(tt2[,sortBy]),] # Sort by t-statistic

    # Select the up/down-regulated gene sets
    mouse.upreg.hits = unique(tt3[dim(tt3)[1]:(dim(tt3)[1] - thresh), "MGI.symbol"])
    mouse.downreg.hits = unique(tt3[1:thresh,"MGI.symbol"])

    # Select the background gene set
    mouse.bg = unique(tt3$MGI.symbol)

    # Do EWCE analysis
    #full_res_up = bootstrap.enrichment.test(sct_data=sct_data,annot=annot,mouse.hits=mouse.upreg.hits,mouse.bg=mouse.bg,reps = reps,annotLevel=1,geneSizeControl=FALSE)
    #full_res_down = bootstrap.enrichment.test(sct_data=sct_data,annot=annot,mouse.hits=mouse.downreg.hits,mouse.bg=mouse.bg,reps = reps,annotLevel=1,geneSizeControl=FALSE)
    full_res_up = EWCE::bootstrap.enrichment.test(sct_data=sct_data,hits=mouse.upreg.hits,bg=mouse.bg,reps = reps,annotLevel=annotLevel,geneSizeControl=FALSE,genelistSpecies=ttSpecies,sctSpecies=sctSpecies)
    full_res_down = EWCE::bootstrap.enrichment.test(sct_data=sct_data,hits=mouse.downreg.hits,bg=mouse.bg,reps = reps,annotLevel=annotLevel,geneSizeControl=FALSE,genelistSpecies=ttSpecies,sctSpecies=sctSpecies)

    joint_results = rbind(cbind(full_res_up$results,Direction="Up"),cbind(full_res_down$results,Direction="Down"))
    hit.cells.up = full_res_up$hit.cells
    hit.cells.down = full_res_down$hit.cells
    bootstrap_data.up = full_res_up$bootstrap_data
    bootstrap_data.down = full_res_down$bootstrap_data

    
    result = list(joint_results=joint_results,hit.cells.up=hit.cells.up,hit.cells.down=hit.cells.down,bootstrap_data.up=bootstrap_data.up,bootstrap_data.down=bootstrap_data.down)
    return(list(joint_results=joint_results,hit.cells.up=hit.cells.up,hit.cells.down=hit.cells.down,bootstrap_data.up=bootstrap_data.up,bootstrap_data.down=bootstrap_data.down))
    
    #write.csv(result,file=paste("Results_ewce_expression/Tables/p120/", 'p120_', gsub('.csv','','Vent_Med_White'), '.csv', sep=''), quote=F,row.names = F)
}
