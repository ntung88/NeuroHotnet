library('igraph')

#Library of functions for implementing the HotNet FDR algorithm

#Runs the algorithm, returns list of vectors where each vector is a connected component detected as significant
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 results are the averages)
#struct: regions x regions structural data (result of diffusion)
#alpha: significance level
#beta: false discovery rate
#k: max component size tested
#delta: threshold for reduced influence graph (struct)
#trace: should it print results for each threshold?
FDR <- function(dats, struct,alpha,beta,k,delta,numtrials,trace=FALSE) {
  
  #Works for when dats isn't a list too
  if(!is.list(dats)) {dats = list(dats)}
  p = nrow(dats[[1]])
  n = length(dats)
  
  #Average results for each matrix in dats
  model = matrix(integer(length = p^2),nrow=p)
  counts = integer(length = p)
  pvals = integer(length = p)
  samps = integer(length = p)
  for(i in 1:n) {
    En = H(dats[[i]],struct,delta)
    samp = SubCounts(En)
    mc = MC(dats[[i]],struct,numtrials,delta,pthresh=samp)
    model = model + En/n
    samps = samps + samp/n
    counts = counts + mc$counts/n
    pvals = pvals + mc$pvals/n
  }
  
  #Construct betas to decrease with component size as specified in hotnet
  betas = integer(k)
  sofar = 0
  for(i in k:2) {
    betas[i] = beta/2^(k-i+1)
    sofar = sofar + betas[i]
  }
  betas[1] = beta - sofar
  
  #Look for smallest S that is statistically significant and satisfies FDR condition
  s = 1
  while(s <= k) {
    if (pvals[s] <= alpha/k && samps[s] >= counts[s]/betas[s]) {
      break
    }
    s = s + 1
  }
  if(s == k+1) {s=Inf}
  
  subs = SubNetworks(model)
  
  if(trace) {
    if(s > k+1) {
      print('No s satisfies conditions!')
    } else {
      print(sprintf('FOUND: s=%d',s))
    }
    print(sprintf('Alpha: %f, Beta: %f, K: %d, Threshold: %f',alpha,beta,k,delta))
    print(data.frame(
      s = 2:k,
      Num_Components = samps[2:k],
      pvals = pvals[2:k],
      Exp_Num_Comp = counts[2:k],
      betas = betas[2:k]
    ))
  }
  
  grs = Filter(function(x) length(x) >= s, subs)
  return(list('groups'=grs,
              'pvals'=formatC(pvals[lengths(grs)],format='e',digits=2)))
}

#Same as FDR but returns binary matrix form of results, not list of groups
#Note this only works for a single observation matrix, doesn't average across multiple
FDR_Mat <- function(dat, struct,alpha,beta,k,delta,numtrials) {
  
  En = H(dat,struct,delta)
  stopifnot(min(En) >= 0 && isSymmetric(En))
  
  betas = integer(k)
  sofar = 0
  for(i in k:2) {
    betas[i] = beta/2^(k-i+1)
    sofar = sofar + betas[i]
  }
  betas[1] = beta - sofar
  
  samp = SubCounts(En)
  mc = MC(dat,struct,numtrials,delta,pthresh=samp)
 
  #Look for smallest S that is statistically significant and satisfies FDR condition
  s = 1
  while(s <= k) {
    if (mc$pvals[s] <= alpha/k && samp[s] >= mc$counts[s]/betas[s]) {
      break
    }
    s = s + 1
  }
  if(s == k+1) {
    print('No s satisfies conditions!')
    return(NaN)
  } else {
    print(sprintf('FOUND: s=%d',s))
  }
  
  #zero out non connected components
  subs = SubNetworks(En)
  for(group in Filter(function(x) length(x) < s, subs)) {
    En[group,] = 0
    En[,group] = 0
  }
  
  return(En>0)
}

#Shuffle rows of dat (observation data), used for null hypothesis
GenNull <- function(dat) {
  return(dat[sample(nrow(dat)),])
}

#Compute enhanced influence model using correlation from observations
#S: Observation data
#GI: structural matrix/influence graph
#delta: threshold to get reduced influence graph
H <- function(S,GI,delta) {
  #CHANGE TO MULTILAYER LATER?
  scalars = cor(t(S))
  En = GI * scalars
  En[which(En < delta)] = 0
  return(En)
}

#Run monte carlo simulation for pvalues and counts of number of connected components of varying sizes
#Returns list of vectors 'counts' and 'pvals'. The value at index [i] of these vectors is the count/pval for components of size i
#S: Observation data
#GI: structural matrix/influence graph
#n: number of mc trials to run
#delta: threshold to get reduced influence graph
#pthresh: vector of values derived from real observation data, used as threshold to calculate pvalues
MC <- function(S,GI,n,delta,pthresh=NA) {
  nregions = nrow(S)
  counts = integer(nregions)
  pvals = integer(nregions)
  for(i in 1:n) {
    nullH = H(GenNull(S),GI,delta)
    subs = SubCounts(nullH)
    counts = counts + subs
    if(!is.na(pthresh)) {
      pvals = pvals + (subs >= pthresh)
    }
  }
  
  if(is.na(pthresh)) {
    return(list('counts' = counts/n,'pvals' = NA))
  }
  return(list('counts' = counts/n,'pvals' = pvals/n))
}

#Utility function to get counts for every possible size of connected component in m
#as above, value at index i is count for components of size i
SubCounts <- function(m,pretty=FALSE) {
  stopifnot(nrow(m)==ncol(m))
  grps = SubNetworks(m)
  counts = integer(nrow(m))
  for(g in grps) {
    counts[1:length(g)] = counts[1:length(g)] + 1
  }
  if(pretty) {
    return(data.frame(
      s = 1:30,
      Num_Components = counts[1:30]
    ))
  }
  return(counts)
}

#Utility function to get connected groups in m
SubNetworks <- function(m) {
  g  <- graph.adjacency(m>0,mode = 'undirected')
  clu <- components(g)
  return(Filter(function(x) length(x) > 1, groups(clu)))
}