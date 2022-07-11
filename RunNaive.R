#Script to run SC-Naive algorithm
library('R.matlab')
source('tasks.R')
library('psych')
source('naive.R')
type = 'Time Courses/'
subjects = list.files(path = type,pattern = '^tmmatrx')

tcs = list()
for(i in 1:length(subjects)) {
  load(sprintf('%s%s',type,subjects[i]))
  rownames(tm.matrx) <- NULL
  tcs[[i]] = tm.matrx
}

threshold = 1e-79
nv = naive(tcs,threshold,mu=0.5)
print(num2name(nv$groups))
