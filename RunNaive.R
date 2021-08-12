setwd('~/Desktop/NeuroHotNetHD')
# source('heat.R')
# source('FDR.R')
library('R.matlab')
source('tasks.R')
type = 'Time Courses/'
subjects = list.files(path = type,pattern = '^tmmatrx')
# subjects = subjects[1:2]

threshold = 1e-180

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

res = naive(tcs,threshold)
print(res$groups)
# print(res$pvals)