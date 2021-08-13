library('igraph')

#Library of functions for implementing the HotNet FDR algorithm

#Runs the algorithm with Higgins, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is the pvalue of each component in 'groups'
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
FDR <- function(dats,struct,alpha,k,numtrials,delta,pr=FALSE,trace=FALSE,mu=NULL,dats2=NULL,beta=NULL) {
  
  #Works for when dats isn't a list too
  if(!is.list(dats)) {dats = list(dats)}
  p = nrow(dats[[1]])
  n = length(dats)
  
  #Average results for each matrix in dats
  #an collect average across all edges for testing of components later
  model = matrix(integer(length = p^2),nrow=p)
  aggcor = c()
  for(i in 1:n) {
    loccor = clean(cor(t(dats[[i]])))
    aggcor = c(aggcor, mean(loccor,na.rm=TRUE))
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
        slice = sampcor[grs[[i]],grs[[i]]]
        data = c(data, mean(slice,na.rm=TRUE))
        if(!is.null(dats2)) {
          sampcor2 = clean(cor(t(dats2[[j]])))
          slice2 = sampcor2[grs[[i]],grs[[i]]]
          data2 = c(data2, mean(slice2,na.rm=TRUE))
        }
      }
      if(!is.null(mu)) {
        res = t.test(data,mu)
      } else if(!is.null(dats2)) {
        res = t.test(data,data2,paired = pr)
      } else {
        res = t.test(data,aggcor,paired = pr)
      }
      aggres[i] = res$p.value
    }
  }
  
  return(list('groups'=grs,
              'pvals'=formatC(aggres,format='e',digits=2)))
}

#Runs the algorithm with Higgins, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is the pvalue of each component in 'groups'
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 average covariance is used)
#struct: regions x regions structural data (result of diffusion)
#alpha: significance level
#k: max component size tested
#tuning: nu parameter in paper, controls overall sparsity
#pr: should a paired (vs unpaired) t.test be run when testing components
#mu: if non-null, the constant mean under null hypothesis that the mean of dats 
#for each component is t-tested against
#dats2: if non-null AND MU IS NULL, the second set of data for two sample t-test
SiGGM <- function(dats,struct,tuning,pr=FALSE,mu=NULL,dats2=NULL) {
  
  #Works for when dats isn't a list too
  if(!is.list(dats)) {dats = list(dats)}
  p = nrow(dats[[1]])
  n = length(dats)
  
  #Average covariance for each matrix in dats
  #an collect average across all edges for testing of components later
  aggcov = matrix(integer(length = p^2),nrow=p)
  aggcor = c()
  for(i in 1:n) {
    y = t(MyScale(dats[[i]]))
    loccov = (1/nrow(y))*t(y)%*%y
    aggcov = aggcov + loccov/n
    loccor = clean(cor(t(dats[[i]])))
    aggcor = c(aggcor, mean(loccor,na.rm = TRUE))
  }
  
  
  hig = SCFC(matrix(integer(length = p*ncol(dats[[1]])),nrow=p),struct,cov_init=aggcov,c0=tuning)
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
        slice = sampcor[grs[[i]],grs[[i]]]
        data = c(data, mean(slice,na.rm=TRUE))
        if(!is.null(dats2)) {
          sampcor2 = clean(cor(t(dats2[[j]])))
          slice2 = sampcor2[grs[[i]],grs[[i]]]
          data2 = c(data2, mean(slice2,na.rm=TRUE))
        }
      }
      if(!is.null(mu)) {
        res = t.test(data,mu)
      } else if(!is.null(dats2)) {
        res = t.test(data,data2,paired = pr)
      } else {
        res = t.test(data,aggcor,paired = pr)
      }
      aggres[i] = res$p.value
    }
  }
  
  return(list('groups'=grs,
              'pvals'=formatC(aggres,format='e',digits=2)))
}

#Utility function to turn Nan into 0 and then make diagonal Nan (so diagonal can be
#ignored when taking mean)
clean<- function(mat) {
  mat[which(is.na(mat))] = 0
  diag(mat) = NaN
  return(mat)
}

# FDR <- function(dats, struct,alpha,beta,k,delta,numtrials,mode=c('info','mat'),trace=FALSE) {
#   
#   #Works for when dats isn't a list too
#   if(!is.list(dats)) {dats = list(dats)}
#   p = nrow(dats[[1]])
#   n = length(dats)
#   
#   #Average results for each matrix in dats
#   model = matrix(integer(length = p^2),nrow=p)
#   counts = integer(length = p)
#   pvals = integer(length = p)
#   samps = integer(length = p)
#   for(i in 1:n) {
#     En = H(dats[[i]],struct,delta)
#     samp = SubCounts(En)
#     mc = MC(dats[[i]],struct,numtrials,delta,pthresh=samp)
#     model = model + En/n
#     samps = samps + samp/n
#     counts = counts + mc$counts/n
#     pvals = pvals + mc$pvals/n
#   }
#   
#   #Construct betas to decrease with component size as specified in hotnet
#   betas = integer(k)
#   sofar = 0
#   for(i in k:2) {
#     betas[i] = beta/2^(k-i+1)
#     sofar = sofar + betas[i]
#   }
#   betas[1] = beta - sofar
#   
#   #Look for smallest S that is statistically significant and satisfies FDR condition
#   s = 1
#   while(s <= k) {
#     if (pvals[s] <= alpha/k && samps[s] >= counts[s]/betas[s]) {
#       break
#     }
#     s = s + 1
#   }
#   if(s == k+1) {s=Inf}
#   
#   subs = SubNetworks(model)
#   
#   if(trace) {
#     if(s > k+1) {
#       print('No s satisfies conditions!')
#     } else {
#       print(sprintf('FOUND: s=%d',s))
#     }
#     print(sprintf('Alpha: %f, Beta: %f, K: %d, Threshold: %f',alpha,beta,k,delta))
#     print(data.frame(
#       s = 2:k,
#       Num_Components = samps[2:k],
#       pvals = pvals[2:k],
#       Exp_Num_Comp = counts[2:k],
#       betas = betas[2:k]
#     ))
#   }
#   
#   grs = Filter(function(x) length(x) >= s, subs)
#   mode = match.arg(mode)
#   if(mode=='info') {
#     return(list('groups'=grs,
#                 'pvals'=formatC(pvals[lengths(grs)],format='e',digits=2)))
#   } else if (mode=='mat') {
#     if(s > k+1) {
#       return(NaN)
#     }
#     for(group in Filter(function(x) length(x) < s, subs)) {
#       model[group,] = 0
#       model[,group] = 0
#     }
#     return(model>0)
#   }
# }

#Same as FDR but returns binary matrix form of results, not list of groups
#Note this only works for a single observation matrix, doesn't average across multiple
# FDR_Mat <- function(dat, struct,alpha,beta,k,delta,numtrials) {
#   
#   En = H(dat,struct,delta)
#   stopifnot(min(En) >= 0 && isSymmetric(En))
#   
#   betas = integer(k)
#   sofar = 0
#   for(i in k:2) {
#     betas[i] = beta/2^(k-i+1)
#     sofar = sofar + betas[i]
#   }
#   betas[1] = beta - sofar
#   
#   samp = SubCounts(En)
#   mc = MC(dat,struct,numtrials,delta,pthresh=samp)
#  
#   #Look for smallest S that is statistically significant and satisfies FDR condition
#   s = 1
#   while(s <= k) {
#     if (mc$pvals[s] <= alpha/k && samp[s] >= mc$counts[s]/betas[s]) {
#       break
#     }
#     s = s + 1
#   }
#   if(s == k+1) {
#     print('No s satisfies conditions!')
#     return(NaN)
#   }
#   
#   print(sprintf('FOUND: s=%d',s))
#   return(En)
#   
#   # #zero out non connected components
#   # subs = SubNetworks(En)
#   # for(group in Filter(function(x) length(x) < s, subs)) {
#   #   En[group,] = 0
#   #   En[,group] = 0
#   # }
#   # 
#   # return(En>0)
# }

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
  En = GI * S
  En[which(abs(En) < delta)] = 0
  return(En)
}

# H <- function(S,GI,delta) {
#   #CHANGE TO MULTILAYER LATER?
#   scalars = cor(t(S))
#   En = GI * scalars
#   En[which(En < delta)] = 0
#   return(En)
# }

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
    # nullS = S[sample(nrow(S)),sample(ncol(S))]
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

# MC <- function(S,GI,n,delta,pthresh=NA) {
#   nregions = nrow(S)
#   counts = integer(nregions)
#   pvals = integer(nregions)
#   for(i in 1:n) {
#     nullH = H(GenNull(S),GI,delta)
#     subs = SubCounts(nullH)
#     counts = counts + subs
#     if(!is.na(pthresh)) {
#       pvals = pvals + (subs >= pthresh)
#     }
#   }
#   
#   if(is.na(pthresh)) {
#     return(list('counts' = counts/n,'pvals' = NA))
#   }
#   return(list('counts' = counts/n,'pvals' = pvals/n))
# }

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