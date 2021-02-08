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
