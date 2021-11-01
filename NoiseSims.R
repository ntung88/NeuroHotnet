setwd('~/Desktop/NeuroHotnetHD/')
library('R.matlab')
source('heat.R')
source('num2name.R')
source('simulate.R')
source('FDR.R')
source('naive.R')

subj = "tmmatrx_102513.rda"
mat = readMat('whole_brain_AAL2.mat')$connectivity
gamma = 30
n=308
m=2

#### Simulation 1, Nathan's approach, using the true time course mean/sd and true structural
Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
En = H(Linv,matrix(1,nrow=120,ncol=120),0.13)
En = zero.out(En,3)
real = SubNetworks(En)
# linvthresh = 0.17

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
for(j in 1:m) {
  print(j)
  
  tcs = list()
  for(i in 1:n) {
    tcs[[i]] = simulate(subj, En,noisevar = 35000)
  }
  
  fdr = FDR(tcs,En,0.1,5,1000,0.00039)
  # nai = naive(tcs,4e-3)
  res = SiGGM(tcs,En,0.000291)
  
  neurrates = c(neurrates,recovery_rate(fdr$groups,real))
  higrates = c(higrates,recovery_rate(res$groups,real))
  # naiverates = c(naiverates,recovery_rate(nai$groups,real))
}

neurrate = mean(neurrates)
higrate = mean(higrates)
naiverate = mean(naiverates)

# if(trace) {
#   print('NeuroHotnet Results')
#   print(neurrate)
#   print('SiGGM Results')
#   print(higrate)
#   print('Naive Results')
#   print(naiverate)
# }
