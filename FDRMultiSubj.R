setwd('~/Desktop/NeuroHotNetHD')
# source('heat.R')
source('FDR.R')
library('R.matlab')
#102614

type = 'Time Courses/'
aggfdr = matrix(integer(length = 120^2),nrow=120)
subjects = list.files(path = type,pattern = '^tmmatrx')
subjects = subjects[1:100]

gamma = 30
dat.dir = "~/Desktop/2020NeuroHotnet/Onsets"
ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)

### Load the structural weights
mat = readMat('whole_brain_AAL2.mat')$connectivity

Linv = heat(mat,gamma,0,weighted=TRUE, trans = FALSE)

for (subj in subjects) {
  print(sprintf('Starting Subject %s',subj))
  load(sprintf('%s%s',type,subj))
  tc=tm.matrx
  rownames(tc) <- NULL
  
  # weighted, normed, no trans
  deltas = rev(seq(from=1e-5,to=1.2e-4,length.out=15))
  # unweighted,  normed, no trans
  # deltas = rev(seq(from=2e-4,to=9e-4,length.out=15))
  # unweighted,  no normed, no trans
  # deltas = rev(seq(from=2e-4,to=5e-4,length.out=15))
  # deltas = rev(seq(from=1e-5,to=1e-4,length.out=15))
  # normed
  # deltas = seq(from=4.5e-3,to=8e-3,length.out=15)
  # deltas = rev(seq(from=1.8e-4,to=2.5e-4,length.out=12))
  for(d in deltas){
    res = FDR(tc,Linv,0.05,0.05,10,d,500,mode='m',trace = TRUE)
    if(!is.nan(res)) {
      aggfdr = aggfdr + res
      break
    }
  }
}
aggfdr = aggfdr/length(subjects)