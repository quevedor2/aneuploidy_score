####################################################
#### cytoband data: ucsc.[GENOMEBUILD].cytoband ####
# converts cytoband data to R data
pkg <- file.path("~", "git", "aneuploidy_score")
cytoband_files <- list.files(file.path(pkg, "data-raw"), pattern = "cytoband")
cytoband_files <- setNames(cytoband_files, 
                           gsub("^ucsc\\.(.*?)\\..*", "\\1", cytoband_files))

readAndFormat <- function(f_name){
  f <- cytoband_files[[f_name]]
  cband <- read.table(file.path(pkg, "data-raw", f), header=FALSE, check.names = FALSE,
                      stringsAsFactors = FALSE, sep="\t")
  colnames(cband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
  return(cband)
}

ucsc.hg19.cytoband <- readAndFormat('hg19')
usethis::use_data(ucsc.hg19.cytoband, overwrite=T)

ucsc.hg38.cytoband <- readAndFormat('hg38')
usethis::use_data(ucsc.hg38.cytoband, overwrite=T)

ucsc.mm9.cytoband <- readAndFormat('mm9')
usethis::use_data(ucsc.mm9.cytoband, overwrite=T)

ucsc.mm10.cytoband <- readAndFormat('mm10')
usethis::use_data(ucsc.mm10.cytoband, overwrite=T)

ucsc.mm39.cytoband <- readAndFormat('mm39')
usethis::use_data(ucsc.mm39.cytoband, overwrite=T)

#####################################
#### Demo seg file (CCLE: 639-V) ####
# CCLE name: 639V_URINARY_TRACT
# SNP6 Filename: VAULT_p_NCLE_DNAAffy14_GenomeWideSNP_6_F05_698718
# Gender: M	
# Site primary: urinary_tract	
# Histology: carcinoma	
# Histology subtype: transitional_cell_carcinoma
# converts cytoband data to R data
seg <- read.table(file=file.path("~", "git", "aneuploidy_score", "data-raw", "demo.seg"), sep="\t", 
                  stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
seg <- GenomicRanges::makeGRangesFromDataFrame(seg, keep.extra.columns = TRUE)
colids <- c('nMajor', 'nMinor', 'nAraw', 'nBraw', 'TCN', 'seg.mean')
seg <- seg[,colids]

usethis::use_data(seg, overwrite=T)

