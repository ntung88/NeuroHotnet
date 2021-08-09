setwd('~/Desktop/NeuroHotNetHD')
# source('heat.R')
source('MyScale.R')
source('SCFCpath.R')
source('SCFC.R')
source('EBIC.R')
library('R.matlab')

type = 'Time Courses/'
agghig = matrix(integer(length = 120^2),nrow=120)
subjects = list.files(path = type,pattern = '^tmmatrx')
subjects = subjects[1:10]

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
  
  tuning = 0.01
  hig = SCFC(MyScale(tc),Linv,c0=tuning)
  agghig = agghig + hig$Covariance[[1]]
  
  # tuning = seq(from=0.004,to=0.009,length.out=10)
  # tuning = seq(from=0.002,to=0.009,length.out=20)
  # tuning = c(seq(from=0.02,to=0.04,length.out=2),seq(from=0.08,to=0.15,length.out=8))
  # 0.004
  # tuning=c(0.0072)
  # tuning = seq(from=0.005,to=0.03,length.out=10)
  # hig = SCFCpath(MyScale(tc),Linv,nulist=tuning)
  # agghig = agghig + hig$Covariance[[1]]
  
  
  # print(SubNetworks(hig$Covariance))
}
agghig = agghig/length(subjects)