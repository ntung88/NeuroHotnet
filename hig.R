setwd('~/Desktop/NeuroHotNetHD')
# source('heat.R')
source('MyScale.R')
source('SCFCpath.R')
# source('FDR.R')
source('EBIC.R')
#102614
load('Time Courses/tmmatrx_102513.rda')
tc = tm.matrx
rownames(tc) <- NULL
gamma = 30
Linv = heat(gamma,0,weighted=FALSE, normed = TRUE, trans = FALSE)
# Linv = Linv/max(Linv)
# Linv[3,] <- 0
# Linv[,3] <- 0

# tuning = 0.012
# hig = SCFC(MyScale(tc),Linv,c0=tuning,a0_init = 100, b0_init = 0.1)

# tuning = seq(from=0.004,to=0.009,length.out=10)
# tuning = seq(from=0.002,to=0.009,length.out=20)
# tuning = c(seq(from=0.02,to=0.04,length.out=2),seq(from=0.08,to=0.15,length.out=8))
# 0.004
# tuning=c(0.0072)
# hig = SCFCpath(MyScale(tc),Linv,nulist=tuning)
# 
# print(SubNetworks(hig$Covariance))

# weighted, normed, no trans
# deltas = rev(seq(from=1e-5,to=1.2e-4,length.out=15))
# unweighted,  normed, no trans
deltas = rev(seq(from=2e-4,to=9e-4,length.out=15))
# unweighted,  no normed, no trans
# deltas = rev(seq(from=2e-4,to=5e-4,length.out=15))
# deltas = rev(seq(from=1e-5,to=1e-4,length.out=15))
# normed
# deltas = seq(from=4.5e-3,to=8e-3,length.out=15)
for(d in deltas){
  res = FDR(tc,Linv,0.05,0.05,20,d,500,trace = TRUE)
  if(length(res$groups)>0) {
    print('FOUND')
    print(res$groups)
    break
  }
}