MyScale <- function(tc) {
  tctil=tc-matrix(rowMeans(tc),nrow=nrow(tc),ncol=ncol(tc))
  dividend = matrix(rowSums(tctil^2)^.5,nrow=nrow(tctil),ncol=ncol(tctil))
  dividend[which(dividend == 0)] <- Inf
  return(tctil/dividend)
}