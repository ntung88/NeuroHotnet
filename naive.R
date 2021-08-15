

naive <- function(dats,delta,pr=FALSE,mu=NULL,dats2=NULL) {
  #collect average correlation for each component across all subjects
  #and test for pvalue
  p = nrow(dats[[1]])
  n = length(dats)
  
  aggcor = c()
  arr = array(0, c(p,p,n))
  arr2 = array(0, c(p,p,n))
  for(i in 1:n) {
    loccor = clean(cor(t(dats[[i]])))
    arr[,,i] = loccor
    if (!is.null(dats2)) {
      arr2[,,i] = clean(cor(t(dats2[[i]])))
    } else {
      aggcor = c(aggcor, mean(loccor,na.rm=TRUE))
    }
  }
  
  aggres = matrix(1,p,p)
  for (i in 1:p) {
    for (j in i:p) {
      if(i == j) {next}
      data = arr[i,j,]
      if(!is.null(mu)) {
        res = t.test(data,mu)
      } else if(!is.null(dats2)) {
        res = t.test(data,arr2[i,j,],paired=pr)
      } else {
        res = t.test(data,aggcor,paired = pr)
      }
      aggres[i,j] = res$p.value
      aggres[j,i] = res$p.value
    }
  }
  
  return(list('groups'=Filter(function(x) length(x) > 1, SubNetworks(aggres<delta)),'pvals' = aggres))
}