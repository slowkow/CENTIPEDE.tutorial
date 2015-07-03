#' A dataset ready to use with CENTIPEDE.
#' 
#' @format A list with two items:
#' \describe{
#'   \item{mat}{A numeric matrix with counts of read start sites.
#'    Each row corresponds to a genomic site flanked by 100 bp. A single row
#'    contains information about the positive strand and the negative strand.
#'    So, for example, this dataset uses a transcription factor motif of length
#'    13 bp and flanks of length 100 bp. Therefore, the matrix has
#'    2 * (100 + 13 + 100) = 426 columns.
#'   }
#'   \item{regions}{A dataframe.
#'    Each row corresponds to a genomic site, in the same order as the
#'    matrix of counts. This has information from a run of FIMO using the 13 bp
#'    motif for STAT4.
#'   }
#' }
"cen"

#' A vector of mean conservation scores for the PWM sites in "cen"
#' 
#' @format A vector that corresponds to the data in "cen". Each item in the
#' vector is the mean conservation score across 100 verterbrates computed
#' with phastCons at UCSC.
#' @source \url{http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/}
"site_cons"