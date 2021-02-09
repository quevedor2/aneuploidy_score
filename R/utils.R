#' listData
#' @description Lists all the datasets included in the package and their
#' category
#' @export
#'
#' @examples
#' listData()
listData <- function(){
  ## Data structure listing all included data files
  datas <- list("cytoband"=list("hg19"="ucsc.hg19.cytoband",
                                "hg38"="ucsc.hg38.cytoband",
                                "mm10"="ucsc.mm10.cytoband",
                                "mm9"="ucsc.mm9.cytoband",
                                "mm39"="ucsc.mm39.cytoband"),
                "test"=list("test"="test"))
  
  ## simple function to formate the print statement for cat()
  datas_print <- sapply(names(datas), function(category){
    items <- datas[[category]]
    items_print <- paste(sapply(names(items), function(item_id){
      item = datas[[category]][[1]]
      return(paste0("\t", item_id, "\t|\t", items[[item_id]], "\n"))
    }), collapse = "")
    
    category_print <- paste0(" > ", category, ":\n", items_print)
    return(category_print)
  })
  
  cat(paste(datas_print, collapse="\n"))
  return(NULL)
}

#' Validates if the seg is compatible with this R package
#'
#' @param seg a Seg file dataframe
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @return LOGICAL
#'
#' @examples
#' #.validateSeg(data.frame("Chrom"=1, "start"=1, "end"=2))
#' #.validateSeg(data.frame("chr"=1, "X"=1, "loc.end"=2))
.validateSeg <- function(seg){
  seg_check <- TRUE
  
  ## GRanges compatibility check
  seg_check <- tryCatch({
    makeGRangesFromDataFrame(seg[1,,drop=FALSE])
    TRUE
  }, error=function(e){ 
    FALSE
  })
  return(seg_check)
}

#' Gets fraction of intervals that make up each GR object
#' @description Takes a GRanges object (ideally and designed for
#' a single chromosome or chromosome arm) and founds out what fraction
#' of the total chr/chr_arm each interval makes up of it. It does this
#' for all intervals with a copy number value (non-NA in cn_col), as 
#' well as with NA's removed (no CN value (NA) in cn_col).
#'
#' @param gr GRanges object of a chromosome or chromosome arm
#' @param cn_col Column in GRanges object containing CN values [Default='CN']
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges width
#' @return 2-column data frame of chr_fractions with NA's included 
#' and NA's excluded
.getChrarmFractions <- function(gr, cn_col='CN'){
  na_incl <- round(width(gr) / sum(width(gr)),3)
  na_excl <- if(any(is.na(mcols(gr)[,cn_col]))){
    round(width(gr) / 
            sum(width(gr[-which(is.na(mcols(gr)[,cn_col])),])),3)
  } else {
    na_incl
  }
  data.frame("na_incl" = na_incl,
             "na_excl" = na_excl)
}