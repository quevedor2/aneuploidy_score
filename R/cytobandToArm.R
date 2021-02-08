#' cytobandToArm: Reduces cytoband data to chromosome arms
#' @description Takes the data([CYTOBAND]) information and splits the 
#' cytobands into their individual p-arm, q-arm or centromere positions.
#' Then it will reduce and clean up the data before returning it as a list
#' of chromosomal positions with 3 rows: 'p', 'q', and 'cen'.
#' @param cytoband cytoband data from data(ucsc.[GENOMEBUILD].cytoband)
#'
#' @return
#' @export
#'
#' @examples
#' data(ucsc.hg19.cytoband)
#' cytobandToArm(ucsc.hg19.cytoband)
cytobandToArm <- function(cytoband){
  # Get a list of the main chromosomes
  chrs <- unique(gsub("_.*", "", cytoband$chrom))
  if(any(grepl("chrM|chrUn", chrs))) chrs <- chrs[-grep("chrUn|chrM", chrs)]
  chrs <- chrs[na.omit(match(paste0("chr", c(1:100, "X", "Y")), chrs))]
  
  # Split cytobands by chromosomes
  cytoband_chr <- split(cytoband, f=cytoband$chrom)

  ## Cycle through each chromosome to reduce to p-, q-arm or centromere
  chr_arms <- lapply(cytoband_chr, function(cyt){
    ## Identify which cytobands are p-, q-arm or centromere
    arms <- as.character(factor(gsub("^([pq]).*", "\\1", cyt$name), c("p", "q", "cen")))
    arms[grep("cen", cyt$gieStain)] <- 'cen'
    levels(arms) <- c('p', 'cen', 'q')
    
    ## Reduce the arm start-end position and cytoband ranges
    cyt_arms <- lapply(split(cyt, f=arms), function(cyt_arm){
      cyt_arm[1,'chromEnd'] <- (cyt_arm[nrow(cyt_arm),'chromEnd']-1)
      cyt_arm[1,'name'] <- paste0(cyt_arm[1,'name'], "-", cyt_arm[nrow(cyt_arm),'name'])
      cyt_arm[1,'gieStain'] <- NA
      return(cyt_arm[1,,drop=FALSE])
    })
    cyt_arms <- as.data.frame(do.call(rbind, cyt_arms))
    
    return(cyt_arms)
  })
  
  ## Reorder chromosomes and return list of p-, q-arm, centromere
  #chr_ord <- order(factor(names(chr_arms), levels=paste0("chr", c(1:22, "X", "Y"))))
  chr_ord <- match(chrs, names(chr_arms))
  chr_arms <- chr_arms[chr_ord]
  return(chr_arms)
}
