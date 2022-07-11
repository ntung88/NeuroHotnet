#Script to run recovery-rate simulation under noise
source('heat.R')
library('Matrix')
source('naive.R')

gamma = 30
n=308
m=5
noisevar=25000

naiparam = 1e-3
resparam = 0.00035
higparam = 0.00025

for (batch in 2) {
  print(sprintf('BATCH: %d',batch))
  load(sprintf('hypo%d.rda',batch))
  Linv = heat(hypo,gamma,0,weighted=TRUE, trans = TRUE)
  real = SubNetworks(hypo)
  
  err <- function( x, y) { sum(lengths(setdiff( union(x, y), intersect(x, y))))}
  
  recovery_rate <- function(estimate,true) {
    correct = 0
    for(gr in true) {
      if (Position(function(x) err(gr,x) <= 2, estimate, nomatch = 0) > 0) {
        correct = correct + 1
      }
    }
    return(correct/length(true))
  }
  
  neurrates = c()
  higrates = c()
  naiverates = c()
  oldhigrates = c()
  for(j in 1:m) {
    print(sprintf('TRIAL: %d',j))
    
    tcs = list()
    for(i in 1:n) {
      tcs[[i]] = simFC(hypo,noisevar)
    }
    
    fdr = NeurHot(tcs,Linv,0.1,10,1000)
    nai = naive(tcs,naiparam)
    res = SiGGM(tcs,Linv,resparam)
    hig = SiGGM(tcs,hypo,higparam)
    
    neurrates = c(neurrates,recovery_rate(fdr$groups,real))
    higrates = c(higrates,recovery_rate(res$groups,real))
    naiverates = c(naiverates,recovery_rate(nai$groups,real))
    oldhigrates = c(oldhigrates,recovery_rate(hig$groups,real))
    
    if(j %% 25 == 0) {
      results = list(naiverates,oldhigrates,higrates,neurrates)
      save(results,file=sprintf('rates_%d_intermediate_2022.rda',batch-1))
    }
  }
  
  results = list(naiverates,oldhigrates,higrates,neurrates)
  save(results,file=sprintf('rates_%d_2022.rda',batch-1))
  
  print(lapply(results, mean))
}
