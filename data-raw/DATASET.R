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
