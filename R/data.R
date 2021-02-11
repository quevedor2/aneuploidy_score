#' hg19 UCSC cytoband table
#'
#' A data table containing the cytoband and centromere information
#' for hg19 genome
#'
#' @format A data frame with 862 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome (e.g. chr1)}
#'   \item{chromStart}{Chromosome cytoband start position (0-based)}
#'   \item{chromEnd}{Chromosome cytoband end position (0-based)}
#'   \item{name}{Name of cytoband}
#'   \item{gieStain}{Giemsa staining}
#'   ...
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/}
"ucsc.hg19.cytoband"

#' hg38 UCSC cytoband table
#'
#' A data table containing the cytoband and centromere information
#' for hg38 genome
#'
#' @format A data frame with 1433 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome (e.g. chr1)}
#'   \item{chromStart}{Chromosome cytoband start position (0-based)}
#'   \item{chromEnd}{Chromosome cytoband end position (0-based)}
#'   \item{name}{Name of cytoband}
#'   \item{gieStain}{Giemsa staining}
#'   ...
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/}
"ucsc.hg38.cytoband"

#' mm10 UCSC cytoband table
#'
#' A data table containing the cytoband and centromere information
#' for mm10 genome
#'
#' @format A data frame with 403 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome (e.g. chr1)}
#'   \item{chromStart}{Chromosome cytoband start position (0-based)}
#'   \item{chromEnd}{Chromosome cytoband end position (0-based)}
#'   \item{name}{Name of cytoband}
#'   \item{gieStain}{Giemsa staining}
#'   ...
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/}
"ucsc.mm10.cytoband"

#' mm39 UCSC cytoband table
#'
#' A data table containing the cytoband and centromere information
#' for mm39 genome
#'
#' @format A data frame with 121 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome (e.g. chr1)}
#'   \item{chromStart}{Chromosome cytoband start position (0-based)}
#'   \item{chromEnd}{Chromosome cytoband end position (0-based)}
#'   \item{name}{Name of cytoband}
#'   \item{gieStain}{Giemsa staining}
#'   ...
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/}
"ucsc.mm39.cytoband"

#' mm9 UCSC cytoband table
#'
#' A data table containing the cytoband and centromere information
#' for mm9 genome
#'
#' @format A data frame with 403 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome (e.g. chr1)}
#'   \item{chromStart}{Chromosome cytoband start position (0-based)}
#'   \item{chromEnd}{Chromosome cytoband end position (0-based)}
#'   \item{name}{Name of cytoband}
#'   \item{gieStain}{Giemsa staining}
#'   ...
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/}
"ucsc.mm9.cytoband"

#' Demo .seg file
#'
#' A .seg file from CCLE's 639-V urinary carcinoma cell line, processed using ASCAT
#' with the allele-specific and total copy numbers represented, as well as the L2R.
#' The data frame was converted into a GenomicRanges object.
#'
#' @format A data frame with 384 ranges and 6 metadata columns variables:
#' \describe{
#'   \item{nMajor}{Haplotype A specific copy number}
#'   \item{nMinor}{Haplotype B specific copy number}
#'   \item{nAraw}{Haplotype A specific copy number (raw)}
#'   \item{nBraw}{Haplotype B specific copy number (raw)}
#'   \item{TCN}{Total copy number}
#'   \item{seg.mean}{Log 2 Ratio (L2R) segment mean}
#'   ...
#' }
#' @source \url{https://portals.broadinstitute.org/ccle}
"seg"
