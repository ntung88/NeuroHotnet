tasks <- function(subj) {
  
  TR = 0.72

  ### Load the time courses
  load(sprintf('Time Courses/%s',subj))
  rownames(tm.matrx) <- NULL
  onset.files = dir(sprintf('Onsets/%sEVs',gsub('.*_(.+)\\..*','\\1',subj)), pattern = "*", full.names = TRUE)
  
  ### Read in the onsets for each of the 5 tasks
  ### Right hand time indices
  rh.ons = round(read.delim(onset.files[5], header = FALSE, sep = "")[,1]/TR)
  rh.ind = c(rh.ons[1]:(rh.ons[1]+16), rh.ons[2]:(rh.ons[2]+16))
  
  ### Left foot time indices
  lf.ons = round(read.delim(onset.files[2], header = FALSE, sep = "")[,1]/TR)
  lf.ind = c(lf.ons[1]:(lf.ons[1]+16), lf.ons[2]:(lf.ons[2]+16))
  
  ### Left hand time indices
  lh.ons = round(read.delim(onset.files[3], header = FALSE, sep = "")[,1]/TR)
  lh.ind = c(lh.ons[1]:(lh.ons[1]+16), lh.ons[2]:(lh.ons[2]+16))
  
  ### Right foot time indices
  rf.ons = round(read.delim(onset.files[4], header = FALSE, sep = "")[,1]/TR)
  rf.ind = c(rf.ons[1]:(rf.ons[1]+16), rf.ons[2]:(rf.ons[2]+16))
  
  ### tongue time indices
  t.ons = round(read.delim(onset.files[7], header = FALSE, sep = "")[,1]/TR)
  t.ind = c(t.ons[1]:(t.ons[1]+16), t.ons[2]:(t.ons[2]+16))
  
  ### tm.matrx is an R by T patrix of time courses
  ### Extract the time course for right hand
  tm.matrx.rh = tm.matrx[,rh.ind]
  
  ### Extract the time course for left hand
  tm.matrx.lh = tm.matrx[,lh.ind]
  
  ### Extract the time course for right foot
  tm.matrx.rf = tm.matrx[,rf.ind]
  
  ### Extract the time course for left foot
  tm.matrx.lf = tm.matrx[,lf.ind]
  
  ### Extract the time course for tongue
  tm.matrx.t = tm.matrx[,t.ind]
  
  ### Compute functional connectivity graphs for each of the two matrices
  ### above to describe task specific functional connectivity
  
  return(list('RHand'=tm.matrx.rh,'LHand'=tm.matrx.lh,'RFoot'=tm.matrx.rf,'LFoot'=tm.matrx.lf,'Tongue'=tm.matrx.t))
}