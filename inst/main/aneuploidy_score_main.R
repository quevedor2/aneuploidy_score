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

## > Centromere and chromosome-size detection
# Uses the hardcoded position from 'Shukla et al.' or the positions from BioConductor
# https://github.com/pascalduijf/CAAs_1/blob/master/scripts/CAAs_from_cell_lines.py


getChrLength <- function(){
  require(BSgenome.Hsapiens.UCSC.hg19)
  chr.lengths = seqlengths(Hsapiens)[1:24]
  chr.len.gr <- makeGRangesFromDataFrame(data.frame("chrom"=names(chr.lengths),
                                                    "loc.start"=rep(1, length(chr.lengths)),
                                                    "loc.end"=chr.lengths))
  chr.len.gr$cum.end <- cumsum(as.numeric(end(chr.len.gr)))
  chr.len.gr$cum.start <- chr.len.gr$cum.end - (end(chr.len.gr) -1)
  chr.len.gr$cum.mid <- (chr.len.gr$cum.start + ((chr.len.gr$cum.end - chr.len.gr$cum.start)/2))
  return(chr.len.gr)
}

##############
#### Main ####
chr.size.dat <- getChrLength()


seg <- read.table(file.path(PDIR, seg_file), sep="\t", header=TRUE, 
                  stringsAsFactors = FALSE, check.names = FALSE)