library(MASS)
library(miscTools)
library(igraph)
library(psych)
library(ggplot2)
library(DescTools)
library(glasso)
library(R.matlab)

#mat: matrix representation of the graph to run diffusion on
#gamma: flow rate parameter
#thresh: pre-threshold before running diffusion
#weighted: consider magnitudes of entries in mat (weights on edges)?
#trans: convert to transition matrix after diffusion (row and col sums 1)?
#plot:plot stages?
heat <- function(mat,gamma,thresh,weighted=TRUE,trans=FALSE,plot=FALSE) {
  
  #For weighted diffusion we don't thresh old
  if (weighted) {
    adj = mat
  } else {
    adj = mat > thresh
  }
  
  #Tools to vizualize the process if plot is set to True
  if(plot) {
    if(weighted) {
      plotadj = final.mat > 2000
    } else {
      plotadj = adj
    }
    diag(plotadj) <- 0
    gr = graph_from_adjacency_matrix(plotadj,mode='undirected')
    V(gr)$name = V(gr)
    isolated = which(degree(gr)==0)
    gr = delete.vertices(gr, isolated)
    pdf('Stage1.pdf')
    plot(gr,vertex.size=8,main=sprintf('Adjacency Matrix Prior to Diffusion (threshold: %d)',thresh))
    dev.off()
  }
  
  #For weighted diffusion we normalize by degree and end up with nonnegative
  #values. For unweighted we just keep the binary adjacency matrix
  
  norm = diag(apply(adj, 2, sum))
  norm = norm^(-1/2)
  norm[which(norm==Inf)] <- 0
  adj = norm %*% adj %*% norm
  
  D = diag(apply(adj, 2, sum))
  A = adj
  L = -A + D
  L = L + gamma*diag(nrow(L))
  
  Linv = solve(L)
  diag(Linv) <- 0
  
  if(trans) {
    sums = rowSums(Linv)
    sums[which(sums==0)] <- Inf
    p = nrow(Linv)
  
    for(i in 1:p) {
      for(j in i:p) {
        avg = (Linv[i,j]/sums[i] + Linv[j,i]/sums[j])/2
        Linv[i,j] = avg
        Linv[j,i] = avg
      }
    }
  }

  if(plot) {
    plotL = Linv
    if (weighted) {
      gr = graph_from_adjacency_matrix(plotL>0.0012,mode='undirected')
    } else {
      gr = graph_from_adjacency_matrix(plotL>0.00085,mode='undirected')
    }
    V(gr)$name = V(gr)
    isolated = which(degree(gr)==0)
    gr = delete.vertices(gr, isolated)
    pdf('Stage2.pdf')
    plot(gr,vertex.size=8,main='Structural Matrix Post Diffusion (only most strongly connected edges)')
    dev.off()
  }
  
  return(Linv)
}

#Like heat but randomized, useful for simulations
rheat <- function(mat,gamma,thresh,weighted=TRUE) {
  
  mat = mat[sample(nrow(mat)),sample(ncol(mat))]
  mat = (mat + t(mat))/2
  
  if (weighted) {
    adj = mat > 0
    norm = diag(apply(adj, 2, sum))
    norm = norm^(-1/2)
    norm[which(norm==Inf)] <- 0
    adj = norm %*% mat %*% norm
  } else {
    adj = mat > thresh
  }
  
  D = diag(apply(adj, 2, sum))
  A = adj
  L = -A + D
  L = L + gamma*diag(nrow(L))
  
  Linv = solve(L)
  diag(Linv) <- 0
  
  return(Linv)
}



