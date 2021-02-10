#' Whole genome doubling checker
#' 
#' @description Takes a .seg file data frame and estimates whether
#' the tumor is whole-genome doubled (WGD). This function requries
#' Total Copy Number and will not work using L2R.
#' @param segf .Seg Data [Data frame]
#' @param tcn_col Column ID for CN value [String] 
#' @param threshold The threshold to call diploid
#' (Diploid = wgd_ploidy +/- threshold) (Default: 0.5) [Numeric]
#' @param wgd_gf Genomic fraction above wgd_ploidy
#' to be considered genome-doubled (Default: 0.5) [Numeric]
#' @param wgd_ploidy Base ploidy (Default: 2) [Integer]
#'
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors mcols
#' @importFrom stats weighted.mean
#' @return A named vector
#'  ploidy: The mean CN weighted by segment sizes
#'  WGD: Boolean indicated whether sample is WGD
#'   
#' @export
checkIfWGD <- function(segf, tcn_col,threshold=0.5,
                       wgd_gf = 0.5, wgd_ploidy=2){
  # tcn_col = 'Modal_Total_CN'
  stopifnot(.validateSeg(segf))
  
  ## Set up GenomicRange Objects of cytoband and seg files
  seg_gr <- makeGRangesFromDataFrame(segf, keep.extra.columns = TRUE)
  seqlevelsStyle(seg_gr) <- 'UCSC'
  
  # Validation checks
  cn_stat <- .classifyCN(cn=round(mcols(seg_gr)[,tcn_col],1), ploidy=2, threshold=threshold)
  width_cn <- data.frame("width", width(seg_gr),
                         "CN"=round(mcols(seg_gr)[,tcn_col],1),
                         "CN_stat"=cn_stat)
  
  # Calculate the fraction of the genome for each loss, neut, gain
  total_genome_size <- sum(width_cn$width)
  cn_genome_width <- sapply(split(width_cn, f=width_cn$CN_stat), function(i) sum(i$width))
  cn_genome_fraction <- cn_genome_width / total_genome_size
  
  # Check if WGD is present
  WGD_status <- as.logical(cn_genome_fraction['GAIN'] > wgd_gf)
  ploidy <- weighted.mean(width_cn$CN, width_cn$width/1000000)
  
  return(c("ploidy"=ploidy, "WGD"=WGD_status))
}