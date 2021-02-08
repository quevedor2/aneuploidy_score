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

## > L2R to call a gain or loss
# "using a threshold of >0.2 for amplification and <−0.2 for deletion"
l2r_threshold_gain <- 0.2
l2r_threshold_loss <- -0.2

## > Remove segments that overlap the centromere
# " For each segment in the SCNA file, if the segment intersected the centromere 
# positions it was discarded"
filter_centromere <- FALSE

## > Centromere and chromosome-size detection
# Uses the hardcoded position from 'Shukla et al.' or the positions from BioConductor
# https://github.com/pascalduijf/CAAs_1/blob/master/scripts/CAAs_from_cell_lines.py

##############
#### Main ####
seg <- read.table(file.path(PDIR, seg_file), sep="\t", header=TRUE, 
                  stringsAsFactors = FALSE, check.names = FALSE)

seg_samples <- split(seg, f=seg$Sample)
segf <- seg_samples[[1]]

library(AneuploidyScore)
data("ucsc.hg19.cytoband")
cytoband <- ucsc.hg19.cytoband
cytoarm <- cytobandToArm(cytoband)


#' Title
#'
#' @param segf 
#' @param cytoband 
#' @importFrom GenomicsRanges makeGRangesFromDataFrame
#' @importFrom GenomicsRanges seqnames
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @return
#' @export
#'
#' @examples
getCAA <- function(segf, cytoband, tcn_col,
                   filter_centromere=FALSE){
  stopifnot(.validateSeg(segf))
  
  ## Set up GenomicRange Objects of cytoband and seg files
  cyto_gr <- lapply(cytoarm, makeGRangesFromDataFrame, keep.extra.columns=TRUE)
  seg_gr <- makeGRangesFromDataFrame(segf, keep.extra.columns = TRUE)
  seqlevelsStyle(seg_gr) <- 'UCSC'
  seg_chr <- split(seg_gr, f=seqnames(seg_gr))
  
  ## Goes through chr by chr, 
  seg_cyto_chr <- lapply(names(seg_chr), function(chr_id){
    segc <- seg_chr[[chr_id]]
    cytoc <- sort(cyto_gr[[chr_id]])
    
    ## Create a GRanges object with all unique intervals between segc and cytoc
    starts <- sort(c(start(segc), start(cytoc)))
    ends <- sort(c(end(segc), end(cytoc)))
    combc <- GRanges(seqnames=chr_id, 
                     IRanges(start=unique(sort(c(starts, ends[-length(ends)]+1))),
                             end=unique(sort(c(ends, starts[-1]-1)))))
    
    if(filter_centromere){
      ## // TODO 
      # Implement segment removal for centromere overlap
      ## TODO \\
    } else {
      # Map chr-arm to intervals
      cyto_comb_ov <- findOverlaps(cytoc, combc)
      mcols(combc)$arm <- names(cytoc[queryHits(cyto_comb_ov),])
      
      # Map seg values to intervals
      segc_comb_ov <- findOverlaps(segc, combc)
      mcols(combc)$CN <- NA
      mcols(combc[subjectHits(segc_comb_ov),])$CN <- mcols(segc[queryHits(segc_comb_ov),])[,tcn_col]
    }
    
    ## Handle intervals that have no CN value (NA; i.e. telomeric ends)
    if(any(is.na(combc$CN))) {
      combc$gf_na <- round(width(combc)/sum(width(combc)),3) ## Genome fraction with NA values
      combc$UID <- paste(combc$arm, combc$CN, sep="_")  ## Sets a unique arm_CN value
      combc_na <- combc[which(is.na(combc$CN)),]
      combc <- combc[-which(is.na(combc$CN)),] ## Removes intervals with no CN values (NA)
      combc$gf <- round(width(combc)/sum(width(combc)),3) ## Genome fraction without NA values
    }
    
    ## Reduce intervals where there is no CN change
    # e.g. c(4,4,4,2,4,4) => list(c(4,4,4), c(2), c(4,4))
    uid_int <- as.integer(factor(combc$UID))
    combc_arms <- split(combc, cumsum(c(TRUE, diff(uid_int) != 0)))
    combc_arms <- as(lapply(combc_arms, function(ca){
      ca_red <- reduce(ca)
      mcols(ca_red) <- mcols(ca)[1,]
      return(ca_red)
    }), "GRangesList")
    
    
    return(unlist(combc_arms))
  })
  names(seg_cyto_chr) <- names(seg_chr)
  
  return(seg_cyto_chr)
}


checkIfWGD <- function(){
  
}

