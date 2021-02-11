#' Reduce Arm/Seg level CN to an aneuploidy score
#'
#' @description Takes the Arm and Seg level CN classification (Loss (-1),
#' Neut (0), Gain(1)) from the getCAA() function and reduces this data 
#' into a single matrix of chromosomal-arm level value (e.g. chr1p: 1)
#' 
#' @param segf_caa Output from getCAA(), GRanges of individual 
#' chromosomes [GRangesList]
#' @param caa_method The chromosomal-arm aneuploidy scoring metric to use: [Vector]
#'  - 'arm': Reduces the arm-level CN classifications (classified from the Weighted Median of the chr-arm CN) to a single value 
#'  - 'seg': Uses a threshold (or weighted median if no threshold given) on the chr-arm CN classifications across all arm-level segs 
#' @param threshold For 'seg' caa_method, how much of a chr-arm must be gained or lost to get -1 or 1 score
#' @param arm_ids [Default c("p", "q")]
#' @param verbose verbose statement
#' 
#' @importFrom stats setNames
#' @importFrom matrixStats weightedMedian
#' @importFrom GenomicRanges width
#' @importFrom S4Vectors expand.grid
#' 
#' @return
#' A matrix of aneuploidy scores [chr-arms x caa_methods]
#' 
#' @export
reduceArms <- function(segf_caa, caa_method=c('arm', 'seg'), threshold=NA, 
                       arm_ids=c('p', 'q'), verbose=FALSE){
  arm_caa_mat <- sapply(caa_method, function(cam){
    if(verbose) cat(paste0("Summarizing aneuploidy score for '", toupper(cam), "'...\n"))
    
    if(toupper(cam) == 'ARM'){
      ## Takes the "Arm-level" CN classes and reports the aneuploidy score in a matrix
      # Arm-level is based off of the 'Weighted Median' CN of each Chr arm
      arms_caa <- lapply(segf_caa, function(s){
        arms <- sapply(split(s$armCNclass, s$arm)[arm_ids], unique)
        suppressWarnings(setNames(as.numeric(as.character(arms)), arm_ids))
      })
      arms_caa <- as.matrix(unlist(arms_caa))
    } else if(toupper(cam) == 'SEG'){
      ## Takes the "Seg-level" CN classes,  and reports the aneuploidy score in a matrix
      # Seg-level reports class for each segment and requires further processing for arm-calls
      arms_caa <- lapply(segf_caa, function(s){
        f_arm <- factor(s$arm, levels=arm_ids)
        
        arms <- sapply(split(s[,'segCNclass'], f_arm), function(sarm){
          if(is.na(threshold)) {
            ## If no threhsold given, just returns the weightedMedian classification
            weightedMedian(sarm$segCNclass, width(sarm))
            
          } else if (threshold >= 0 & threshold <= 1) {
            ## Looks for the CN segments make up a certain size of the chr-arm
            arm_sizes <- sapply(split(sarm, sarm$segCNclass), function(cn_size) {
              sum(width(cn_size)) / sum(width(sarm))
            })
            if(any(arm_sizes > threshold)){
              as.integer(names(arm_sizes[arm_sizes>threshold]))
            } else if(length(arm_sizes) == 0) { 
              NA
            } else {
              0
            }
            
          } else {
            stop("Threshold must be NA or a numeric value from 0-1")
          }
        })
        suppressWarnings(setNames(as.numeric(as.character(arms)), arm_ids))
      })
      arms_caa <- as.matrix(unlist(arms_caa))
    }
    
    return(arms_caa)
  })
  rownames(arm_caa_mat) <- apply(expand.grid(arm_ids, names(segf_caa))[,c(2,1)], 1, 
                                 paste, collapse="_")
  return(arm_caa_mat)
}
