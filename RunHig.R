setwd('~/Desktop/NeuroHotNetHD')
source('heat.R')
# source('FDR.R')
library('R.matlab')
source('tasks.R')
source('MyScale.R')
source('SCFC.R')
type = 'Time Courses/'
subjects = list.files(path = type,pattern = '^tmmatrx')

# nu for all task data
nu = 0.011
# nu = 0.007

#nu for task-separated data
# nu = 0.09

gamma = 30
dat.dir = "~/Desktop/2020NeuroHotnet/Onsets"
ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)

tcs = list()
rh = list()
lh = list()
rf = list()
lf = list()
ton = list()
for(i in 1:length(subjects)) {
  load(sprintf('%s%s',type,subjects[i]))
  rownames(tm.matrx) <- NULL
  tcs[[i]] = tm.matrx
  tasked = tasks(subjects[i])
  rh[[i]] = tasked$RHand
  lh[[i]] = tasked$LHand
  rf[[i]] = tasked$RFoot
  lf[[i]] = tasked$LFoot
  ton[[i]] = tasked$Tongue
}

### Load the structural weights
mat = readMat('whole_brain_AAL2.mat')$connectivity

Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
# sink('HigginsResults.txt')
res = SiGGM(tcs,Linv,nu)
print(res$groups)
# print(res$pvals)
# sink()