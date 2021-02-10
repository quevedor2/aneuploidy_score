#devtools::install_github("quevedor2/aneuploidy_score", ref='master')
devtools::install_github("quevedor2/aneuploidy_score", ref='debug')

library(AneuploidyScore)
## For Calculating Aneuploidy Score
## Similar to the way done by Shukla et al. 2020
## Counts the number of chromosomal arm gains/losses adjusted for ploidy
PDIR <- "/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA/input"
seg_file <- "TCGA_mastercalls.abs_segtabs.fixed.txt"

####################
#### PARAMETERS ####
## > Fraction of the genome (gf) needed to be duplicated (copies) to be called whole-genome doubled
# "if more than half of the autosomal genome had two or more copies of the more 
# frequent (maternal or paternal) allele"
wgd_threshold <- c("gf"=0.5, "copies"=2)

## > Length to call an arm-level aberration
# "Chromosome arm-level SCNAs/aneuploidies (CAAs) were called by scoring 
# individual chromosome arms as gained or lost if ≥0.9 of the arm was gained or lost."
arm_threshold <- 0.9

## > CN reporting method
# Whether a total CN method is being used or a L2R
cn_method <- 'TCN' # TCN or L2R

## > L2R (or TCN) to call a gain or loss
# L2R: "using a threshold of >0.2 for amplification and <−0.2 for deletion"
# TCN: 0.5 based on round() function from
# https://github.com/broadinstitute/Aneuploidy_dependencies/blob/master/make_CCLE_arm_calls.R
threshold <- 0.5

## > Remove segments that overlap the centromere
# " For each segment in the SCNA file, if the segment intersected the centromere 
# positions it was discarded"
filter_centromere <- FALSE

## > Centromere and chromosome-size detection
# Uses the hardcoded position from 'Shukla et al.' or the positions from BioConductor
# https://github.com/pascalduijf/CAAs_1/blob/master/scripts/CAAs_from_cell_lines.py

## > CN value column
# In the seg file, what is the column header to be used for gain/loss/neutral
tcn_col <- 'Modal_Total_CN'

##############
#### Main ####
seg <- read.table(file.path(PDIR, seg_file), sep="\t", header=TRUE, 
                  stringsAsFactors = FALSE, check.names = FALSE, nrows = 1000)

seg_samples <- split(seg, f=seg$Sample)
segf <- seg_samples[[2]]




library(AneuploidyScore)
data("ucsc.hg19.cytoband")
cytoband <- ucsc.hg19.cytoband
cytoarm <- cytobandToArm(cytoband)

wgd_ploidy <- checkIfWGD(segf, tcn_col = tcn_col)
segf_caa <- getCAA(segf, cytoarm, tcn_col=tcn_col, classifyCN=TRUE,
                   ploidy=wgd_ploidy['ploidy'], threshold=threshold)


#' @importFrom matrixStats weightedMean
#' #' @importFrom matrixStats weightedMedian

## Rescale CN by ploidy if using Total Copy Number
if(grepl('tcn', cn_method, ignore.case = TRUE)){
  segf_gr <- unlist(segf_caa)
  
  ## Assign base ploidy to closest to a multi of base2 (e.g. 2,4,6,8)
  ploidy_multi <- cut(wgd_ploidy['ploidy'], breaks=seq(-1,11,by=2), right=FALSE)
  levels(ploidy_multi) <- seq(0,10, by=2)
  
  ploidy_val <- switch(cen_method,
                       "wmean"=weightedMean(segf_gr$CN, (width(segf_gr)/1*10^6)),
                       "wmedian"=weightedMedian(segf_gr$CN, (width(segf_gr)/1*10^6)),
                       "mean"=mean(segf_gr$CN),
                       "median"=median(segf_gr$CN),
                       "multi_base2"=as.integer(as.character(ploidy_multi)),
                       NA)
  
  ## get CN-classes per CN-segment
  segf_gr$deltaCN <- .classifyCN(cn=segf_gr$CN, ploidy=ploidy_val, threshold=0.5)
  
  ## get CN-classes per weighted-median chr-arm
  chrs_segf_gr <- as(lapply(split(segf_gr, f=seqnames(segf_gr)), function(chr_segf){
    arms_segf <- as(lapply(split(chr_segf, f=chr_segf$arm), function(arm_segf){
      arm_segf$armCN <- rep(round(weightedMedian(arm_segf$CN),2), length(arm_segf))
      return(arm_segf)
    }), "GRangesList")
    return(unlist(arms_segf))
  }), "GRangesList")
  segf_gr <- sort(unlist(chrs_segf_gr))
  segf_gr$deltaCNarm <- .classifyCN(cn=segf_gr$armCN, ploidy=ploidy_val, threshold=0.5)
  
  return(segf_gr)
}




