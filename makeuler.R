#Script to make Euler diagram of groups from different algorithms
#sources necessary to define variables neur, res, nv with results
source('RunNeur.R')
source('RunNaive.R')
source('RunsiG.R')
library('Rcpp')
library('eulerr')
library('grid')

fdrgrs = Filter(function(x) length(x) > 3, neur$groups)
higgrs = Filter(function(x) length(x) > 3, res$groups)
nvgrs = Filter(function(x) length(x) > 3, nv$groups)

combos = list()
for(i in 1:length(fdrgrs)) {
  combos[paste('H',i,sep='')] = fdrgrs[i]
}
for(i in 1:length(higgrs)) {
  combos[paste('S',i,sep='')] = higgrs[i]
}
for(i in 1:length(nvgrs)) {
  combos[paste('N',i,sep='')] = nvgrs[i]
}
#hardcode regions with perfect overlap
combos[c('H2','H3','H4','H5','S2','S3','S4','S5')] <- NULL
combos['H2,S2'] = fdrgrs[2]
combos['H3,S3'] = fdrgrs[3]
combos['H4,S4'] = fdrgrs[4]
combos['H5,S5'] = fdrgrs[5]
eul = euler(combos,shape='ellipse')
pal = rainbow(4,alpha=0.7)
#harcode ncom with number of regions with perfect overlap
ncom=4
p = plot(eul,main='Neuro-Hotnet vs SiGGM vs SC Naive',labels=list(cex=0.7),quantities=list(cex=0.7),fills=c(rep(pal[4],length(fdrgrs)-ncom),rep(pal[3],length(higgrs)-ncom),rep(pal[1],length(nvgrs)),rep('#4080FFB3',ncom)),cex=0.1)
lg = legendGrob(labels=c('Neuro-Hotnet','SiGGM','SC Naive'),gp=gpar(col=c(pal[4],pal[3],pal[1]), fill="gray",cex=0.7),pch=c(19,19,19))
g <- gTree(children = gList(p))
pdf('Euler.pdf')
gridExtra::grid.arrange(g, lg, nrow = 2, heights = c(5,1))
dev.off()



