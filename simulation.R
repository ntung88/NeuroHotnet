subj = "tmmatrx_102513.rda"
mat = readMat('whole_brain_AAL2.mat')$connectivity
gamma = 30
n=308

#### Simulation 1, Nathan's approach, using the true time course mean/sd and true structural
# Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
# tcs = list()
# for(i in 1:n) {
#   tcs[[i]] = simulate(subj, Linv)
# }
# d=0.04
# print('Simulation 1 Simple Results')
# fdr = FDR(tcs,Linv,0.05,10,1000,d,trace=TRUE)
# print(fdr$groups)
# print(fdr$pvals)

#### Simulation 2, the data are generated from hypothetical 120 regions. 
#### Regions 1:10 are interconnected with correlation 0.8, 
#### Regions 21:50 are interconnected with correlation 0.4.
#### No effect of structural

# Sigma = matrix(0, 120, 120)
# Sigma[1:10, 1:10] = 0.8
# Sigma[21:50, 21:50] = 0.4
# diag(Sigma) = 1
# image(Sigma)
# 
# tcs = list()
# for(i in 1:n) {
#   #I know this simulation should be independent of any real subject data but
#   #it's only using the subject for row means and standard deviations
#   #to keep scale consistent (so we can keep similar parameters especially for higgins), 
#   #the simulated data is only used for its correlation anyway
#   #so there isn't any real dependence on the subject. The correlation will
#   #be roughly sigma with some noise as intended
#   tcs[[i]] = simulate(subj, Sigma)
# }
# 
# d=0.02
# print('Simulation 2 Simple Results')
# 
# #Can choose either truly no effect of structural or uniform random effect
# #of structural, results are robust and as expected either way!
# # fdr = FDR(tcs,matrix(1,120,120),0.05,10,1000,d,trace=TRUE)
# fdr = FDR(tcs,matrix(runif(120^2),120,120),0.05,10,1000,d,trace=TRUE)
# 
# print(fdr$groups)
# print(fdr$pvals)


# 
# ### Simulation 3, the data are generated from hypothetical 120 regions. 
# #### Regions 1:10 are interconnected with correlation 0.8, 
# #### Regions 21:50 are interconnected with correlation 0.4.
# #### True structural is used as a prior
#Not sure how this one should work, I don't think I was doing it right, leaving commented out for now
# 
# Sigma = matrix(0, 120, 120)
# Sigma[1:10, 1:10] = 0.8
# Sigma[21:50, 21:50] = 0.4
# diag(Sigma) = rep(1, 120)
# Linv = heat(mat,gamma,0,weighted=TRUE, trans = TRUE)
# test_cor = Linv/max(abs(Linv))
# Sigma = Sigma * test_cor
# image(Sigma)
# 
# tcs = list()
# for(i in 1:n) {
#   tcs[[i]] = simulate(subj, Sigma)
# }
# 
# d=0.02
# print('Simulation 3 Simple Results')
# 
# # fdr = FDR(tcs,matrix(1,120,120),0.05,10,1000,d,trace=TRUE)
# fdr = FDR(tcs,matrix(runif(120^2),120,120),0.05,10,1000,d,trace=TRUE)
# 
# print(fdr$groups)
# print(fdr$pvals)
# FDR(dats = sim.dat, struct = diag(rep(1,120)), alpha = 0.05, beta = 0.2,
#     k = 30, numtrials = 30, delta = 0, trace=FALSE)
# 
# #### Further things to try, have more regions highly connected. 
# #### Fewer regions highly connected 
# #### More or less of noise
# #### Time dependence
# 

### Simulation 4, the data are generated from hypothetical 120 regions. 
# #### Regions 1:10 are interconnected with correlation 0.8, 
# #### Regions 21:50 are interconnected with correlation 0.4.
# #### True structural used in algorithm
