#devtools::install_github("quevedor2/aneuploidy_score", ref='master')
#devtools::install_github("quevedor2/aneuploidy_score", ref='debug')
devtools::install_github("quevedor2/aneuploidy_score", ref='dev')

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




library("AneuploidyScore")
data("seg")
data("ucsc.hg19.cytoband")

cytoband <- ucsc.hg19.cytoband
cytoarm <- cytobandToArm(cytoband)

## Cohen-Sharir TCN approach:
tcn_col <- 'TCN'
wgd_ploidy <- checkIfWGD(seg, tcn_col = tcn_col, threshold = 0.5,
                         wgd_gf = 0.5, ploidy_method = 'wmean')

threshold <- 0.5
seg_caa <- getCAA(seg, cytoarm, tcn_col=tcn_col, classifyCN=TRUE,
                   ploidy=wgd_ploidy['ploidy'], threshold=threshold)

arm_threshold <- 0.9
caa <- reduceArms(seg_caa, caa_method=c('arm', 'seg'), threshold=arm_threshold, 
                  arm_ids = c('p', 'q'))
colSums(abs(caa), na.rm = TRUE)

## Shukla's L2R approach:
l2r_col <- 'seg.mean'
wgd_ploidy <- checkIfWGD(seg, tcn_col = l2r_col, threshold = 0.2,
                         wgd_gf = 0.5, ploidy_method = 'input', input_ploidy=0)
wgd_ploidy

threshold <- 0.2
seg_caa <- getCAA(seg, cytoarm, tcn_col=l2r_col, classifyCN=TRUE,
                  ploidy=wgd_ploidy['ploidy'], threshold=threshold, filter_centromere = TRUE)

arm_threshold <- 0.9
caa <- reduceArms(seg_caa, caa_method=c('arm', 'seg'), threshold=arm_threshold, 
                  arm_ids = c('p', 'q'))
colSums(abs(caa), na.rm = TRUE)
