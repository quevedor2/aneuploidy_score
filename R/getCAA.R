#' Seg to CN Chromosomal Fractions
#' 
#' @description Takes a .seg data-frame and finds out for each CN
#' interval, what fraction of the chromosomal arm or chromosome does it
#' occupy. It calculates this metric across all covered regions in (NA-
#' excluded) or the entire chromosomal-arm/chromosome (NA-included). It 
#' also supports the Shukla et al. method of removing segments that 
#' overlap the centromere.
#' 
#' @param segf .Seg Data [Data frame]
#' @param tcn_col Column ID for CN value [String] 
#' @param cytoarm Output of cytobandToArm() function [List of Data frames]
#' @param filter_centromere Whether to include or exclude segments
#' @param classifyCN classify arm-CN into Loss/Neut/Gain [Boolean]
#' @param ploidy base-ploidy for Loss/Neut/Gain estimation [Default=0]
#' @param threshold Threshold around ploidy (+/- threshold) [Default=0.2]
#' @param ... Parameters for .classifyCN() function (ploidy, threshold, verbose)
#' 
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom methods as
#' @return A GRangesList split by chromosome with each
#' interval representing a reduced chromosomal aberration
#' and its associated chromosomal_arm and chromosome fraction
#' with and without NAs factored in
#' 
#' @export
getCAA <- function(segf, cytoarm, tcn_col,
                   filter_centromere=FALSE, classifyCN=FALSE,
                   ploidy=0, threshold=0.2, ...){
  # tcn_col = 'Modal_Total_CN'
  stopifnot(.validateSeg(segf))
  
  ## Set up GenomicRange Objects of cytoband and seg files
  cyto_gr <- lapply(cytoarm, GenomicRanges::makeGRangesFromDataFrame, 
                    keep.extra.columns=TRUE)
  seg_gr <- GenomicRanges::makeGRangesFromDataFrame(segf, keep.extra.columns = TRUE)
  seqlevelsStyle(seg_gr) <- 'UCSC'
  seg_chr <- split(seg_gr, f=seqnames(seg_gr))
  
  ## Goes through chr by chr, 
  seg_cyto_chr <- lapply(names(seg_chr), function(chr_id){
    segc <- seg_chr[[chr_id]]
    cytoc <- sort(cyto_gr[[chr_id]])
    
    ## Create a GRanges object with all unique intervals between segc and cytoc
    starts <- sort(c(GenomicRanges::start(segc), GenomicRanges::start(cytoc)))
    ends <- sort(c(GenomicRanges::end(segc), GenomicRanges::end(cytoc)))
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
    
    ## Get chromosomal arm or whole-chromosome fractions for each CN interval
    .assembleFrac <- function(combc, assemble_method='arm', ...){
      if(assemble_method == 'arm'){
        X_fract <- lapply(split(combc, f=combc$arm), .getChrarmFractions, ...)
        ord <- order(factor(names(X_fract), levels=c("p", "cen", "q")))
        X_fract <- as.data.frame(do.call(rbind, X_fract[ord]))
        colnames(X_fract)[1:2] <- c("CAA_frac_NA", "CAA_frac_nonNA")
      } else {
        X_fract <- .getChrarmFractions(combc, ...)
        colnames(X_fract)[1:2] <- c("Chr_frac_NA", "Chr_frac_nonNA")
      }
      mcols(combc) <- cbind(mcols(combc), X_fract)
      return(combc)
    }
    
    combc <- .assembleFrac(combc, assemble_method='chr') # Chromosome fractions
    combc <- .assembleFrac(combc, assemble_method='arm', classifyCN=classifyCN, 
                           ploidy=ploidy, threshold=threshold) # Chromosome arm fractions
    
    ## Handle intervals that have no CN value (NA; i.e. telomeric ends)
    if(any(is.na(combc$CN))) {
      combc$UID <- paste(combc$arm, round(combc$CN,2), sep="_")  ## Sets a unique arm_CN value
      combc_na <- combc[which(is.na(combc$CN)),]
      combc <- combc[-which(is.na(combc$CN)),] ## Removes intervals with no CN values (NA)
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
    
    combc_arms <- unlist(combc_arms)
    mcols(combc_arms) <- mcols(combc_arms)[,-grep("^UID$", colnames(mcols(combc_arms)))]
    if(classifyCN){
      combc_arms$segCNclass <- .classifyCN(cn = combc_arms$CN, ploidy = ploidy, 
                                           threshold=threshold)
    }
    
    
    return(combc_arms)
  })
  names(seg_cyto_chr) <- names(seg_chr)
  
  return(as(seg_cyto_chr, "GRangesList"))
}

