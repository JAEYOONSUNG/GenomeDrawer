#' gc_ratio #
#'
#' @param sequence A nucleotide sequence.
#' @return A list containing summary information of the processed data.
#' @export
#'

gc_ratio <- function(x) {
  if (!is.character(x) || length(x) > 1)
    stop("single string expected")
  tmp <- tolower(seqinr::s2c(x))
  nA <- sum(tmp == "a")
  nT <- sum(tmp == "t")
  nG <- sum(tmp == "g")
  nC <- sum(tmp == "c")
  if (nA + nT + nG + nC == 0)
    return(NA)
  return((nC + nG)/(nA + nT + nG + nC))
}
