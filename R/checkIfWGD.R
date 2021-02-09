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
  loss_threshold <- wgd_ploidy - threshold
  gain_threshold <- wgd_ploidy + threshold
  if(loss_threshold < 0) stop("Loss threshold is below 0. Whole-genome duplication requires a total copy number (TCN) input, not L2R")
  if(abs(loss_threshold-2) > 2) warning("A Total Copy Number (TCN) method is selected, but the loss threshold is too low to be informative")
  if(abs(gain_threshold-2) > 2) warning("A Total Copy Number (TCN) method is selected, but the gain threshold is very distant from copy-neutral (2)")
  
  # Break down the CN genome into loss/neutral/gain
  width_cn <- data.frame("width"=width(seg_gr), 
                         "CN"=round(mcols(seg_gr)[,tcn_col],1))
  cn_cuts <- cut(width_cn$CN, breaks = c(-0.1, loss_threshold, gain_threshold, 100))
  levels(cn_cuts) <- c("LOSS", "NEUT", "GAIN")
  width_cn$CN_stat <- as.character(cn_cuts)
  
  # Calculate the fraction of the genome for each loss, neut, gain
  total_genome_size <- sum(width_cn$width)
  cn_genome_width <- sapply(split(width_cn, f=width_cn$CN_stat), function(i) sum(i$width))
  cn_genome_fraction <- cn_genome_width / total_genome_size
  
  # Check if WGD is present
  WGD_status <- as.logical(cn_genome_fraction['GAIN'] > wgd_gf)
  ploidy <- weighted.mean(width_cn$CN, width_cn$width/1000000)
  
  return(c("ploidy"=ploidy, "WGD"=WGD_status))
}