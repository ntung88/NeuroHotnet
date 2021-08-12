library('Matrix')

#Function to simulate observation data with row means and standard deviations
#matching the observations for subject 'subj' and 'Linv' as a correlation prior (if not a proper correlation
#matrix it's converted to the closest one)

#subj: string of the form 'tmmatrx_######.rda' representing the subject
#Linv: correlation prior
#exact: boolean determining if simulated data sample correlation exactly matches prior, or if there is noise
simulate <- function(subj,Linv,exact=FALSE) {
  load(sprintf('Time Courses/%s',subj))
  tc=tm.matrx
  rownames(tc) <- NULL
  stds = apply(t(tc),2,sd)
  test_cor = Linv/max(abs(Linv))
  diag(test_cor) = 1
  test_cor = round(test_cor,10)
  # stopifnot(is.positive.definite(test_cor))
  test_cov = nearPD(stds %*% t(stds) * test_cor)$mat
  return(t(mvrnorm(n=ncol(tc),mu=colMeans(t(tc)),Sigma=test_cov,empirical = exact)))
}

