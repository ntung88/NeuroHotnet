# setwd('~/Desktop/NeuroHotnetHD/')

mult = 50

start = Sys.time()

n=308
rs = seq(from=200,to=500,by=mult)
hottimes = vector(mode = "list", length = length(rs))
sigtimes = vector(mode = "list", length = length(rs))

for (i in 1:length(rs)) {
  r = rs[i]
  sc = simSC(0.3,r)
  
  tcs = list()
  for(j in 1:n) {
    tcs[[j]] = simFC(sc,noisevar = 500)
  }
  
  hotnet = microbenchmark(FDR(tcs,sc,0.05,10,1000,0.1),times=10)
  hottimes[[i]] = hotnet$time
  siggm = microbenchmark(SiGGM(tcs,sc,0.11),times=10)
  sigtimes[[i]]= siggm$time
  
  save(hottimes,file='hottimes.rda')
  save(sigtimes,file='sigtimes.rda')
}

# pdf(file='runtimes.pdf')
# plot(rs,hottimes)
# par(new=TRUE)
# plot(rs,sigtimes)
# dev.off()

print(Sys.time()-start)
