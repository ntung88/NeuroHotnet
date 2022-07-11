#Runs SC-naive algorithm, returns list where 'groups' is a list of vectors where each vector is a connected component 
#detected as significant. 'pvals' is a matrix of pvalues for each edge. 'mat' contains average correlations across subjects
#for each edge. 'meancor' gives the average over whole graph
#dats: list of regions x timepoints observation data matrices (or single, if more than 1 average correlation is used)
#epsilon: threshold to induce sparsity and determine partitioning of discovered components
naive <- function(dats,epsilon,pr=FALSE,mu=NULL) {
  #collect average correlation for each component across all subjects
  #and test for pvalue
  p = nrow(dats[[1]])
  n = length(dats)
  
  #collect z-score average of whole graph over all subjects
  aggcor = c()
  arr = array(0, c(p,p,n))
  for(i in 1:n) {
    loccor = fisherz(clean(cor(t(dats[[i]]))))
    arr[,,i] = loccor
    aggcor = c(aggcor, mean(loccor,na.rm=TRUE))
  }
  
  #t-test edge weights against average or specified value mu
  aggres = matrix(1,p,p)
  for (i in 1:p) {
    for (j in i:p) {
      if(i == j) {next}
      data = arr[i,j,]
      if(!is.null(mu)) {
        res = t.test(data,mu=fisherz(mu),alternative = 'greater')
      } else {
        res = t.test(data,aggcor,paired = pr)
      }
      aggres[i,j] = res$p.value
      aggres[j,i] = res$p.value
    }
  }
  
  return(list('groups'=Filter(function(x) length(x) > 1, SubNetworks(aggres<epsilon)),'pvals' = aggres,
              'mat' = rowMeans(arr,dims = 2),'meancor'=rowMeans(arr,dims=2)))
}