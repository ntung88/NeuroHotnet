library(MASS)
library(miscTools)
library(igraph)
library(psych)
library(ggplot2)
library(DescTools)
library(glasso)
library(R.matlab)

### Load the task onset files
### for more information on onsets and tasks see 
### https://protocols.humanconnectome.org/HCP/3T/task-fMRI-protocol-details.html
### the data we are using are based on the motor task
setwd('~/Desktop/NeuroHotnetHD')

#Function to load connectivity matrix and run weighted diffusion process on it
heat <- function(gamma,thresh,weighted=TRUE,normed=TRUE,trans=FALSE,plot=FALSE) {
  dat.dir = "~/Desktop/2020NeuroHotnet/Onsets"
  ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)
  
  ### Load the structural weights 
  final.mat = readMat('whole_brain_AAL2.mat')$connectivity
  
  #For weighted diffusion we don't threshold
  if (weighted) {
    # adj = final.mat > 0
    adj = final.mat
  } else {
    adj = final.mat > thresh
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
  if(normed) {
    norm = diag(apply(adj, 2, sum))
    norm = norm^(-1/2)
    norm[which(norm==Inf)] <- 0
    adj = norm %*% final.mat %*% norm
  }
  
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
rheat <- function(gamma,thresh,weighted=TRUE) {
  dat.dir = "~/Desktop/2020NeuroHotnet/Onsets"
  ons.file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)
  
  ### Load the structural weights 
  final.mat = readMat('whole_brain_AAL2.mat')$connectivity
  #Shuffle and resymmetrize
  final.mat = final.mat[sample(nrow(final.mat)),sample(ncol(final.mat))]
  final.mat = (final.mat + t(final.mat))/2
  
  if (weighted) {
    adj = final.mat > 0
    norm = diag(apply(adj, 2, sum))
    norm = norm^(-1/2)
    norm[which(norm==Inf)] <- 0
    adj = norm %*% final.mat %*% norm
  } else {
    adj = final.mat > thresh
  }
  
  D = diag(apply(adj, 2, sum))
  A = adj
  L = -A + D
  L = L + gamma*diag(nrow(L))
  
  Linv = solve(L)
  diag(Linv) <- 0
  
  return(Linv)
}



