setwd('~/Desktop/NeuroHotnet/')
library('R.matlab')
source('heat.R')
source('num2name.R')
source('simulate.R')
source('algos.R')
source('naive.R')
library('microbenchmark')

#Script to evaluate how runtime scales with number of regions of interest

# Generate uniform random SC matrix with specified dimensions
# density: fraction of matrix entries that will be nonzero
# roi: number of regions of interest
simSC <- function(density,roi) {
  vec = runif(floor(density * choose(roi,2)))
  vec = c(vec,integer(choose(roi,2)-length(vec)))
  stopifnot(length(vec) == choose(roi,2))
  vec = sample(vec)
  return(symMatrix(vec,roi))
}

# Sample fMRI time course data from a prior and specified noise level
# density: fraction of matrix entries that will be nonzero
# roi: number of regions of interest
simFC <- function(prior,noisevar) {
  load('Time Courses/tmmatrx_102513.rda')
  tc=tm.matrx
  rownames(tc) <- NULL
  stds = apply(t(tc),2,sd)
  stds = rep(stds[1:10],nrow(prior)/10)
  test_cor = prior/max(abs(prior))
  diag(test_cor) = 1
  test_cor = round(test_cor,10)
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

mult = 50

start = Sys.time()

n=308
rs = seq(from=200,to=500,by=mult)
hottimes = vector(mode = "list", length = length(rs))
sigtimes = vector(mode = "list", length = length(rs))

for (i in 1:length(rs)) {
  r = rs[i]
  sc = simSC(0.3,r)
  
  tcs = list()
  for(j in 1:n) {
    tcs[[j]] = simFC(sc,noisevar = 500)
  }
  
  hotnet = microbenchmark(FDR(tcs,sc,0.05,10,1000,0.1),times=10)
  hottimes[[i]] = hotnet$time
  siggm = microbenchmark(SiGGM(tcs,sc,0.11),times=10)
  sigtimes[[i]]= siggm$time
  
  save(hottimes,file='hottimes.rda')
  save(sigtimes,file='sigtimes.rda')
}

# pdf(file='runtimes.pdf')
# plot(rs,hottimes)
# par(new=TRUE)
# plot(rs,sigtimes)
# dev.off()

print(Sys.time()-start)