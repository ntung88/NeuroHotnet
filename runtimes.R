setwd('~/Desktop/NeuroHotnetHD/')
library('R.matlab')
source('heat.R')
source('num2name.R')
source('simulate.R')
source('FDR.R')
source('naive.R')
library('microbenchmark')

simSC <- function(density,roi) {
  vec = runif(floor(density * choose(roi,2)))
  vec = c(vec,integer(choose(roi,2)-length(vec)))
  stopifnot(length(vec) == choose(roi,2))
  vec = sample(vec)
  return(symMatrix(vec,roi))
}

simFC <- function(prior,noisevar) {
  load('Time Courses/tmmatrx_102513.rda')
  tc=tm.matrx
  rownames(tc) <- NULL
  stds = apply(t(tc),2,sd)
  stds = rep(stds[1:10],nrow(prior)/10)
  test_cor = prior/max(abs(prior))
  diag(test_cor) = 1
  test_cor = round(test_cor,10)
  # stopifnot(is.positive.definite(test_cor))
  test_cov = nearPD(stds %*% t(stds) * test_cor)$mat
  means = colMeans(t(tc))
  means = rep(means[1:10],ncol(prior)/10)
  dat = t(mvrnorm(n=ncol(tc),mu=means,Sigma=test_cov))
  if (!is.null(noisevar)) {
    noise = t(mvrnorm(n=ncol(tc),mu=rep(0,nrow(prior)),Sigma=diag(noisevar,nrow=nrow(prior))))
    dat <- dat + noise
  }
  return(dat)
}

n=308
r = 20

sc = simSC(0.3,r)

tcs = list()
for(i in 1:n) {
  tcs[[i]] = simFC(sc,noisevar = 500)
}

hotnet = microbenchmark(FDR(tcs,sc,0.05,10,1000,0.1),times=1)
siggm = microbenchmark(SiGGM(tcs,sc,0.11),times=1)

print(hotnet)
print(siggm)
