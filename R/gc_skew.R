#' gc_skew #
#'
#' @param sequence A nucleotide sequence.
#' @return A list containing summary information of the processed data.
#' @export
#'

gc_skew <- function(x) {
  if (!is.character(x) || length(x) > 1)
    stop("single string expected")
  tmp <- tolower(seqinr::s2c(x))
  nC <- sum(tmp == "c")
  nG <- sum(tmp == "g")
  if (nC + nG == 0)
    return(NA)
  return((nC - nG)/(nC + nG))
}
