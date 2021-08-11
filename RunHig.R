setwd('~/Desktop/NeuroHotNetHD')
# source('heat.R')
# source('FDR.R')
library('R.matlab')
#102614
type = 'Time Courses/'
aggfdr = matrix(integer(length = 120^2),nrow=120)
subjects = list.files(path = type,pattern = '^tmmatrx')

nu = 0.011

gamma = 30
dat.dir = "~/Desktop/2020NeuroHotnet/Onsets"
ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)

tcs = list()
for(i in 1:length(subjects)) {
  load(sprintf('%s%s',type,subjects[i]))
  rownames(tm.matrx) <- NULL
  tcs[[i]] = tm.matrx
}

### Load the structural weights
mat = readMat('whole_brain_AAL2.mat')$connectivity

Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
# sink('HigginsResults.txt')
res = SiGGM(tcs,Linv,0.05,10,nu)
print(res$groups)
print(res$pvals)
# sink()