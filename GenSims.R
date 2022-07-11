#Simulate roi x roi hypothetical structural matrix with specific density
#entries generated from Uniform(0,1)
simSC <- function(density,roi) {
  vec = runif(floor(density * choose(roi,2)))
  vec = c(vec,integer(choose(roi,2)-length(vec)))
  stopifnot(length(vec) == choose(roi,2))
  vec = sample(vec)
  return(symMatrix(vec,roi))
}

simNet <- function(sizes,roi) {
  samples = lapply(sizes,function (x) runif(x))
}

#Simulate hypothetical fMRI observations with roi same as prior
#row means equivalent to subject 102513, and specific variance of 
#multivariate gaussian noise
simFC <- function(prior,noisevar) {
  #Hardcode mult to be greatest divisor of nrows(prior) that is <= 120
  mult = 120
  load('Time Courses/tmmatrx_102513.rda')
  tc=tm.matrx
  rownames(tc) <- NULL
  stds = apply(t(tc),2,sd)
  stds = rep(stds[1:mult],nrow(prior)/mult)
  test_cor = prior/max(abs(prior))
  diag(test_cor) = 1
  test_cor = round(test_cor,10)
  test_cov = nearPD(stds %*% t(stds) * test_cor)$mat
  means = colMeans(t(tc))
  means = rep(means[1:mult],ncol(prior)/mult)
  dat = t(mvrnorm(n=ncol(tc),mu=means,Sigma=test_cov))
  if (!is.null(noisevar)) {
    noise = t(mvrnorm(n=ncol(tc),mu=rep(0,nrow(prior)),Sigma=diag(noisevar,nrow=nrow(prior))))
    dat <- dat + noise
  }
  return(dat)
}