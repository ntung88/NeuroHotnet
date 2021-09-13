setwd('~/Desktop/NeuroHotNetHD')
source('heat.R')
# source('FDR.R')
library('R.matlab')
source('tasks.R')
type = 'Time Courses/'
subjects = list.files(path = type,pattern = '^tmmatrx')
# subjects = subjects[1:2]

threshold = 0.095
# threshold = 0.068

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
# sink('SimpleResults.txt')
res = FDR(tcs,Linv,0.05,10,1000,threshold,trace = TRUE,mu=1.5)
print(res$groups)
print(res$pvals)
# print(num2name(res$groups))
# sink()
# weighted, normed, no trans
# deltas = rev(seq(from=1e-5,to=1.2e-4,length.out=15))
# deltas = rev(seq(from=8.1e-5,to=8.9e-5,length.out=10))
# unweighted,  normed, no trans
# deltas = rev(seq(from=2e-4,to=9e-4,length.out=15))
# unweighted,  no normed, no trans
# deltas = rev(seq(from=2e-4,to=5e-4,length.out=15))
# deltas = rev(seq(from=1e-5,to=1e-4,length.out=15))
# normed
# deltas = seq(from=3e-2,to=1e-1,length.out=15)
# for(d in deltas){
#   res = FDR(tcs,Linv,0.05,10,1000,d,trace = TRUE)
#   if(length(res$groups)>0) {
#     print('FOUND')
#     print(res$groups)
#     print(res$pvals)
#   }
# }