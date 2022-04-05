source('heat.R')
# source('FDR.R')
library('R.matlab')
source('tasks.R')
source('MyScale.R')
source('SCFC.R')

#Script to run siGGM with Difusion on Human Connectome Project data

#Replace with your file structure, subjects should be list of filenames
#you want to load
#ons.file.dir should be list of onset files
setwd('~/Desktop/NeuroHotNetHD')
type = 'Time Courses/'
subjects = list.files(path = type,pattern = '^tmmatrx')
dat.dir = "~/Desktop/2020NeuroHotnet/Onsets"
ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)

nu = 0.011
gamma = 30

tcs = list()
for(i in 1:length(subjects)) {
  load(sprintf('%s%s',type,subjects[i]))
  rownames(tm.matrx) <- NULL
  tcs[[i]] = tm.matrx
}

### Load the structural weights
mat = readMat('whole_brain_AAL2.mat')$connectivity

Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
res = SiGGM(tcs,Linv,nu,mu=0.5)
print(res$groups)
print(res$pvals)