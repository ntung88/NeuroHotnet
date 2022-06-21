setwd('~/Desktop/NeuroHotnet')
source('heat.R')
# source('algos.R')
library('R.matlab')
source('tasks.R')

gamma = 30

#Replace with your file structure, subjects should be list of filenames
#you want to load
#ons.file.dir should be list of onset files
type = 'Time Courses/'
subjects = list.files(path = type,pattern = '^tmmatrx')
# dat.dir = "~/Desktop/NeuroHotNet/Onsets"
# ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)

tcs = list()
for(i in 1:length(subjects)) {
  load(sprintf('%s%s',type,subjects[i]))
  rownames(tm.matrx) <- NULL
  tcs[[i]] = tm.matrx
}

### Load the structural weights
mat = readMat('whole_brain_AAL2.mat')$connectivity
Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
neur = NeurHot(tcs,Linv,0.05,10,1000,7,trace=TRUE)
print(neur$groups)
print(neur$pvals)