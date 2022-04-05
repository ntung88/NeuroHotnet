library('igraph')
source('MyScale.R')
source('SCFC.R')
library('psych')

#Library of functions for Neuro-Hotnet and siGGM with Diffusion

#Runs the NeuroHotnet algorithm, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is the pvalue of each component in 'groups'. 'mat' contains FC estimates for all edges. 's'
#is the found s^*
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 average correlation is used)
#struct: regions x regions structural data (result of diffusion)
#alpha: significance level
#k: max component size tested
#delta: threshold for reduced influence graph (struct)
#numtrials: number of trials for monte carlo
#pr: should a paired (vs unpaired) t.test be run when testing components
#trace: should it print results for each threshold?
#mu: if non-null, the constant mean under null hypothesis that the mean of dats 
#for each component is t-tested against
#dats2: if non-null AND MU IS NULL, the second set of data for two sample t-test
#beta: if non-null the FDR controlling condition is turned on when testing for s^*
NeurHot <- function(dats,struct,alpha,k,numtrials,delta,pr=FALSE,trace=FALSE,mu=0,dats2=NULL,beta=NULL) {
  
  #Works for when dats isn't a list too
  if(!is.list(dats)) {dats = list(dats)}
  p = nrow(dats[[1]])
  n = length(dats)
  
  #Average results for each matrix in dats
  #an collect average across all edges for testing of components later
  model = matrix(integer(length = p^2),nrow=p)
  for(i in 1:n) {
    loccor = clean(cor(t(dats[[i]])))
    diag(loccor) = 0
    model = model + loccor/n
  }
  
  #Run mc for pvals and expected counts
  En = H(model,struct,delta)
  samp = SubCounts(En)
  mc = MC(model,struct,numtrials,delta,pthresh=samp)
  counts = mc$counts
  pvals = mc$pvals
  
  # Construct betas to decrease with component size as specified in hotnet
  if(!is.null(beta)) {
    betas = integer(k)
    sofar = 0
    for(i in k:2) {
      betas[i] = beta/2^(k-i+1)
      sofar = sofar + betas[i]
    }
    betas[1] = beta - sofar
  }
  
  #Look for smallest S that is statistically significant and satisfies FDR condition
  s = 1
  while(s <= k) {
    if (pvals[s] <= alpha/k && (is.null(beta) || samp[s] >= counts[s]/betas[s])) {
      break
    }
    s = s + 1
  }
  if(s == k+1) {s=Inf}
  
  #Filter components of size at least s
  subs = SubNetworks(En)
  grs = Filter(function(x) length(x) >= s, subs)
  
  if(trace) {
    if(s > k+1) {
      print('No s satisfies conditions!')
    } else {
      print(sprintf('FOUND: s=%d',s))
    }
    print(sprintf('Alpha: %f, K: %d, Threshold: %f',alpha,k,delta))
    print(data.frame(
      s = 2:k,
      Num_Components = samp[2:k],
      pvals = pvals[2:k],
      Exp_Num_Comp = counts[2:k]
    ))
  }
  
  #collect average correlation for each component across all subjects
  #and test for pvalue
  aggres = integer(length(grs))
  if(length(grs) > 0) {
    for (i in 1:length(grs)) {
      data = c()
      data2 = c()
      for(j in 1:n) {
        sampcor = clean(cor(t(dats[[j]])))
        slice = abs(fisherz(sampcor[grs[[i]],grs[[i]]]))
        data = c(data, mean(slice,na.rm=TRUE))
        if(!is.null(dats2)) {
          sampcor2 = clean(cor(t(dats2[[j]])))
          slice2 = abs(fisherz(sampcor2[grs[[i]],grs[[i]]]))
          data2 = c(data2, mean(slice2,na.rm=TRUE))
        }
      }
      if(!is.null(dats2)) {
        res = t.test(data,data2,paired = pr)
      } else {
        res = t.test(data,mu=fisherz(mu),paired = pr,alternative = 'greater')
      }
      aggres[i] = res$p.value
    }
  }
  
  return(list('groups'=grs,
              'pvals'=formatC(aggres,format='e',digits=2),
              'mat' = En,
              's' = s))
}

#Runs siGGM with Diffusion, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is the pvalue of each component in 'groups'. 'mat' contains FC estimates for all edges. 's'
#is the found s^*
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 average covariance is used)
#struct: regions x regions structural data (result of diffusion)
#tuning: nu parameter in Higgins paper, controls overall sparsity
#pr: should a paired (vs unpaired) t.test be run when testing components
#mu: if non-null, the constant mean under null hypothesis that the mean of dats 
#for each component is t-tested against
#dats2: if non-null AND MU IS NULL, the second set of data for two sample t-test
#naive: if TRUE turns off effect of structural, SC-naive
SiGGM <- function(dats,struct,tuning,pr=FALSE,mu=0,dats2=NULL,naive=FALSE) {
  
  #Works for when dats isn't a list too
  if(!is.list(dats)) {dats = list(dats)}
  p = nrow(dats[[1]])
  n = length(dats)
  
  #Average covariance for each matrix in dats
  #an collect average across all edges for testing of components later
  aggcov = matrix(integer(length = p^2),nrow=p)
  for(i in 1:n) {
    y = t(MyScale(dats[[i]]))
    loccov = (1/nrow(y))*t(y)%*%y
    aggcov = aggcov + loccov/n
    loccor = clean(cor(t(dats[[i]])))
  }
  
  
  hig = SCFC(matrix(integer(length = p*ncol(dats[[1]])),nrow=p),struct,cov_init=aggcov,c0=tuning,etaInd = !naive)
  En = hig$Covariance
  
  #extract connected components
  subs = SubNetworks(En)
  grs = Filter(function(x) length(x) > 1, subs)

  #collect average correlation for each component across all subjects
  #and test for pvalue
  aggres = integer(length(grs))
  if(length(grs) > 0) {
    for (i in 1:length(grs)) {
      data = c()
      data2 = c()
      for(j in 1:n) {
        sampcor = clean(cor(t(dats[[j]])))
        slice = abs(fisherz(sampcor[grs[[i]],grs[[i]]]))
        data = c(data, mean(slice,na.rm=TRUE))
        if(!is.null(dats2)) {
          sampcor2 = clean(cor(t(dats2[[j]])))
          slice2 = abs(fisherz(sampcor2[grs[[i]],grs[[i]]]))
          data2 = c(data2, mean(slice2,na.rm=TRUE))
        }
      }
      if(!is.null(dats2)) {
        res = t.test(data,data2,paired = pr)
      } else {
        res = t.test(data,mu=fisherz(mu),paired = pr,alternative = 'greater')
      }
      aggres[i] = res$p.value
    }
  }
  
  return(list('groups'=grs,
              'pvals'=formatC(aggres,format='e',digits=2),
              'mat' = En))
}

#zeros out all edges that are not part of a component of size greater than s
# mat: weighted adjacency matrix (matrix form of graph)
# s: minimum component size
# returns mat with zeroed out edges
zero.out <- function(mat,s) {
  bad = Filter(function(x) length(x) <= s, SubNetworks(mat))
  for (net in bad) {
    mat[net,] <- 0
    mat[,net] <- 0
  }
  return(mat)
}

#calculate relative weighted densities of all subnetworks
# mat: weighted adjacency matrix (matrix form of graph)
# returns vector of same length as number of subnetworks. Each element
# is the weighted density normalized by edges and then by density of the whole
#graph
densities <- function(mat) {
  grs = SubNetworks(mat)
  diag(mat) <- 0
  n = length(grs)
  dense = integer(n)
  for (i in 1:n) {
    gr = grs[[i]]
    s = length(gr)
    dense[i] = sum(mat[gr,gr])/2/choose(s,2)
  }
  return(dense / (sum(mat)/2/choose(nrow(mat),2)))
}

#Turns Nan entries into 0's and then makes diagonal Nan (so diagonal can be
#ignored when taking mean)
# mat: numeric matrix
clean<- function(mat) {
  mat[which(is.na(mat))] = 0
  diag(mat) = NaN
  return(mat)
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
  En = GI * S
  En[which(abs(En) < delta)] = 0
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
    nullS = matrix(sample(S),nrow=nrow(S))
    nullH = H(nullS,GI,delta)
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

#Find counts for every possible size of connected component in m
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

#Find connected components in m
SubNetworks <- function(m) {
  g  <- graph.adjacency(m>0,mode = 'undirected')
  clu <- components(g)
  return(Filter(function(x) length(x) > 1, groups(clu)))
}