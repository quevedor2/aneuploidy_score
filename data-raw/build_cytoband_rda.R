# converts cytoband data to R data
pkg <- file.path("~", "git", "aneuploidy_score")
cytoband_files <- list.files(file.path(pkg, "data-raw"), pattern = "cytoband")
cytoband_files <- setNames(cytoband_files, 
                           gsub("^ucsc\\.(.*?)\\..*", "\\1", cytoband_files))

# f_name <- names(cytoband_files)[1]
cbands <- lapply(names(cytoband_files), function(f_name){
  f <- cytoband_files[[f_name]]
  cband <- read.table(file.path(pkg, "data-raw", f), header=FALSE, check.names = FALSE,
                      stringsAsFactors = FALSE, sep="\t")
  colnames(cband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
  print(paste0("Saving: ", file.path(pkg, "data", gsub(".txt$", ".rda", cytoband_files[f_name]))))
  save(cband, file=file.path(pkg, "data", gsub(".txt$", ".rda", cytoband_files[f_name])))
})
