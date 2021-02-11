#' Whole genome doubling checker
#' 
#' @description Takes a .seg file data frame and estimates whether
#' the tumor is whole-genome doubled (WGD). This function requries
#' Total Copy Number and will not work using L2R.
#'
#' @param segf .Seg Data [Data frame]
#' @param tcn_col Column ID for CN value [String] 
#' @param threshold The threshold to call a gain or loss relative
#' to base ploidy (ploidy +/- threshold) (Default: 0.5) [Numeric]
#' @param wgd_gf Genomic fraction above base_ploidy + threshold to be
#' considered genome doubled (Default: 0.5) [Numeric]
#' @param ploidy_method Calculates the base ploidy [Character]
#'  - 'wmean': Segment size weighted CN mean
#'  - 'wmedian': Segment size weighted CN median
#'  - 'mean': CN mean
#'  - 'median': CN median
#'  - 'multi_base2': Nearest base ploidy 2 (e.g. CN=3 -> ploidy of 4)
#' 
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors mcols
#' @importFrom matrixStats weightedMedian
#' @importFrom matrixStats weightedMean
#' @return A named vector
#'  ploidy: Ploidy as estimated by ploidy_method 
#'  WGD: Boolean indicated whether sample is WGD
#'   
#' @export
checkIfWGD <- function(segf, tcn_col,threshold=0.5,
                       wgd_gf = 0.5, ploidy_method="wmean"){
  # tcn_col = 'Modal_Total_CN'
  stopifnot(.validateSeg(segf))
  
  ## Set up GenomicRange Objects of cytoband and seg files
  segf_gr <- makeGRangesFromDataFrame(segf, keep.extra.columns = TRUE)
  seqlevelsStyle(segf_gr) <- 'UCSC'
  segf_gr$CN <- mcols(segf_gr)[,tcn_col]
  
  ## Assign base ploidy to closest to a multi of base2 (e.g. 2,4,6,8)
  wmean_ploidy <- weightedMean(segf_gr$CN, (width(segf_gr)/1*10^6))
  ploidy_multi <- cut(wmean_ploidy, breaks=seq(-1,11,by=2), right=FALSE)
  levels(ploidy_multi) <- seq(0,10, by=2)
  
  ploidy_val <- switch(ploidy_method,
                       "wmean"=wmean_ploidy,
                       "wmedian"=weightedMedian(segf_gr$CN, (width(segf_gr)/1*10^6)),
                       "mean"=mean(segf_gr$CN),
                       "median"=median(segf_gr$CN),
                       "multi_base2"=as.integer(as.character(ploidy_multi)),
                       NA)

  # Calculate the fraction of the genome for each loss, neut, gain
  total_genome_size <- sum(width(segf_gr))
  cn_stat <- .classifyCN(cn=round(segf_gr$CN,1), ploidy=ploidy_val, threshold=threshold)
  cn_genome_width <- sapply(split(width(segf_gr), f=cn_stat), sum)
  cn_genome_fraction <- round(cn_genome_width / total_genome_size,3)
  
  # Check if WGD is present
  WGD_status <- as.logical(cn_genome_fraction['1'] > wgd_gf)

  return(c("ploidy"=ploidy_val, "WGD"=WGD_status))
}