library('igraph')
source('MyScale.R')
source('SCFC.R')
library('psych')

#Library of functions for Neuro-Hotnet and siGGM with Diffusion

#Runs the NeuroHotnet algorithm, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is the pvalue of each component in 'groups'. 'mat' contains FC estimates for all edges. 's'
#is the found s^*. 'totalpval' is the sum of 'pvals'. 'densities' is the average of edge weights in the component normalized
#by the average of all edge weights in the graph
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 average correlation is used)
#struct: regions x regions structural data (result of diffusion)
#alpha: significance level
#k: max component size tested
#numtrials: number of trials for all monte carlo simulations
NeurHot <- function(dats,struct,alpha,k,numtrials) {
  
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
  
  #Make grid of deltas to test
  runs = numtrials
  range = integer(2)
  for(i in 1:runs){
    En = matrix(sample(model),nrow=nrow(model)) * struct
    percentiles = quantile(En, probs = seq(0, 1, 0.0005))
    range = range + c(percentiles[1997],percentiles[2000])
  }
  # for real data use percentiles[1997] to [2000] and 10 numvals if you want to cut out first local min ;)
  range = range/runs
  print(range)
  numvals = 10
  grid = seq(from=range[1], to=range[2], length.out=numvals)
  
  #Run mc for pvals and expected counts
  Ens = lapply(grid,function(z) H(model,struct,z))
  samps = lapply(Ens,function(z) SubCounts(z))
  mc = MC(model,struct,numtrials,grid,pthreshes=samps)
  counts = mc$counts
  pvals = mc$pvals
  
  #Look for smallest S that is statistically significant and satisfies FDR condition
  ss = lapply(pvals,function(z) sstar(z,k,alpha))
  
  #Filter components of size at least s
  subs = lapply(Ens,function(z) SubNetworks(z))
  grs = vector(mode='list',length=numvals)
  for(i in 1:numvals) {
    grs[[i]] = Filter(function(x) length(x) >= ss[i], subs[[i]])
  }
  
  #collect average correlation for each component across all subjects
  #and test for pvalue
  data = lapply(grs, function(z) lapply(z,function(a) c()))
  mus = lapply(grs, function(z) lapply(z,function(a) c()))
  for(j in 1:n) {
    sampcor = clean(cor(t(dats[[j]])))
    nulldat = dats[[j]][sample(nrow(dats[[j]])),]
    nullcor = clean(cor(t(nulldat)))
    for(l in 1:numvals) {
      deltagrs = grs[[l]]
      if(length(deltagrs) > 0) {
        for(i in 1:length(deltagrs)) {
          nullslice = abs(fisherz(nullcor[deltagrs[[i]],deltagrs[[i]]]))
          mus[[l]][[i]] = c(mus[[l]][[i]], mean(nullslice,na.rm=TRUE))
          
          slice = abs(fisherz(sampcor[deltagrs[[i]],deltagrs[[i]]]))
          data[[l]][[i]] = c(data[[l]][[i]], mean(slice,na.rm=TRUE))
        }
      }
    }
  }
  pvals = mus
  for(l in 1:numvals) {
    if(length(data[[l]]) > 0) {
      for(i in 1:length(data[[l]])) {
        pvals[[l]][[i]] = t.test(data[[l]][[i]],mus[[l]][[i]],paired = TRUE,alternative = 'greater')$p.value
      }
    }
  }
  totalps = unlist(lapply(pvals,function(z) sum(unlist(z))))
  totalps[which(totalps == 0)] = Inf
  ind = which(totalps == min(totalps))
  if(length(ind) != 1) {
    print('ISSUE WITH GRID, MORE THAN ONE OR NO MIN')
    print(grid)
  }
  print(sprintf('delta: %f',grid[ind]))
  ind = ind[1]
  
  return(list('groups'=grs[[ind]],
              'pvals'=pvals[[ind]],
              'totalpval' = totalps[[ind]],
              'mat' = Ens[[ind]],
              'densities' = densities(zero.out(Ens[[ind]],ss[[ind]]-1)),
              's' = ss[ind]))
}

#Finds s^* given pvals from monte carlo simulation
#Returns smallest s such that corresponding pvalue is significant using
#Bonferroni correction with alpha and k
sstar <- function(pvals,k,alpha) {
  s = 2
  while(s <= k) {
    if (pvals[s] <= alpha/k) {
      return(s)
    }
    s = s + 1
  }
  return(Inf)
}

#Runs siGGM with Diffusion, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is the pvalue of each component in 'groups'. 'mat' contains FC estimates for all edges.
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 average covariance is used)
#struct: regions x regions structural data (result of diffusion)
#tuning: nu parameter in Higgins paper, controls overall sparsity
#naive: if TRUE turns off effect of structural, SC-naive
SiGGM <- function(dats,struct,tuning,naive=FALSE) {
  
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
  data = lapply(grs, function(z) c())
  mus = lapply(grs, function(z) c())
  for(j in 1:n) {
    sampcor = clean(cor(t(dats[[j]])))
    nulldat = dats[[j]][sample(nrow(dats[[j]])),]
    nullcor = clean(cor(t(nulldat)))
    for(i in 1:length(grs)) {
      nullslice = abs(fisherz(nullcor[grs[[i]],grs[[i]]]))
      mus[[i]] = c(mus[[i]], mean(nullslice,na.rm=TRUE))

      slice = abs(fisherz(sampcor[grs[[i]],grs[[i]]]))
      data[[i]] = c(data[[i]], mean(slice,na.rm=TRUE))
    }
  }
  pvals = mus
  for(i in 1:length(data)) {
    pvals[i] = t.test(data[[i]],mus[[i]],paired = TRUE,alternative = 'greater')$p.value
  }
  
  return(list('groups'=grs,
              'pvals'=pvals,
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
#S: Correlation matrix from observation data
#GI: structural matrix/influence graph
#n: number of mc trials to run
#deltas: list of thresholds to test
#pthreshes: list of vectors of values derived from real observation data, one for each delta, 
#used as threshold to calculate pvalues
MC <- function(S,GI,n,deltas,pthreshes) {
  nregions = nrow(S)
  nvals = length(deltas)
  counts = lapply(vector(mode='list',length=nvals),function(z) integer(nregions))
  pvals = lapply(vector(mode='list',length=nvals),function(z) integer(nregions))
  for(i in 1:n) {
    nullS = matrix(sample(S),nrow=nrow(S))
    for(j in 1:nvals) {
      nullH = H(nullS,GI,deltas[j])
      subs = SubCounts(nullH)
      counts[[j]] = counts[[j]] + subs
      pvals[[j]] = pvals[[j]] + (subs >= pthreshes[[j]])
    }
  }
  counts = lapply(counts,function(x) x/n)
  pvals = lapply(pvals,function(x) x/n)
  return(list('counts' = counts,'pvals' = pvals))
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
  counts[1] = Inf
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