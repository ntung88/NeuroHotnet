#' Fit the model parameters
#'
#' @description  Estimates a sparse inverse covariance matrix where edge-specific shrinkage parameters are informed by the anatomical connectivity information
#'
#' @details Estimates the anatomically-informed functional brain network via an edge-specific lasso penalty. The inverse covariance matrix is estimated via the graphical lasso (Friedman et al., 2007) or a quadratic approximation to the multivariate normal likelihood plus penalty (Hsieh et al., 2011)  The algorithm also estimates the functional network when the anatomical information is ignored.
#'
#'
#'
#'
#' @param y           p (regions of interest) by T (timepoints) matrix-valued timecourse data.
#' @param P           p by p matrix of structural connectivity weights
#' @param mu_init     Mean of the prior distribution on non-anatomical source of variation in functional connectivity
#' @param a0_init     Scale parameter of the gamma prior prior parameter for eta parameter. Ensure a0_init>1
#' @param b0_init     Shape parameter of the gamma prior prior parameter for eta parameter.
#' @param sigmu       Variance parameter for prior on non-anatomical source of variation in functional connectivity
#' @param siglam      Variance parameter for shrinkage parameters.
#' @param maxits      Maximum number of iterations.
#' @param method      Method for updating the inverse covariance matrix. The two options are 'glasso', which implements the graphical lasso of Friedman et al (2007), or 'QUIC' which implements the QUIC
#' @param etaInd      Argument to indicate the inclusion (etaInd=1) or exclusion (etaInd=0) of the structural connectivity information in the estimation procedure.
#' @param eps         Convergence criteria.  Default is 1e-4
#' @param outerits    The number of outer iterations allowed for each update of the inverse covariance matrix.  Value is ignored if method='glasso'
#' @param c0          Tuning parameter controlling overall network sparsity
#' @param cov_init    If an initial p by p covariance matrix is known, specify it here.  The default is NULL.
#' @param alpha_init  If edge-specific shrinkage values are known, please specify here.  This object must be a p by p matrix but is not required for the function to produce estimates.
#'
#'
#' @return
#' A list with components
#' \item{Omega}{Estimated inverse covariance matrix (p by p matrix)}
#' \item{Covariance}{Estimated covariance matrix (p by p matrix)}
#' \item{lambda}{Estimates of the anatomically informed shrinkage factor at each edge (p by p matrix)}
#' \item{Eta}{Scalar valued estimate of the average effect of structural connectivity on functional connectivity}
#' \item{Method}{Method to estimate the inverse covariance matrix}
#' \item{iters}{Number of iterations until convergence criteria reached}
#' \item{Mu}{Estimate of the non-anatomical source of variation in functional connectivity at each edge}
#' \item{LogLike}{Penalized Log-likelihood at convergence}
#' \item{del}{The change in the objective function at covergence.}
#' \item{BIC}{The Bayesian Information Criterion.}
#'
#' @references Higgins, Ixavier A., Suprateek Kundu, and Ying Guo. Integrative Bayesian Analysis of Brain Functional Networks Incorporating Anatomical Knowledge. arXiv preprint arXiv:1803.00513 (2018).
#' @references Jerome Friedman, Trevor Hastie and Robert Tibshirani (2007). Sparse inverse covariance estimation with the lasso. Biostatistics 2007. http://www-stat.stanford.edu/~tibs/ftp/graph.pdf
#' @references Cho-Jui Hsieh, Matyas A. Sustik, Inderjit S. Dhillon, Pradeep Ravikumar. Sparse Inverse Covariance Matrix Estimation Using Quadratic Approximation. Advances in Neural Information Processing Systems, vol. 24, 2011, p. 2330–2338.
#' @examples
#'   fit<-SmallWorld(10,.15,100)
#'   y=as.matrix(fit$Data)
#'   covdat=cov(y)
#'   omegdat=solve(covdat)
#'   locs=which(abs(omegdat)>quantile(omegdat,probs=.8))
#'   temp=matrix(runif(100,0,1),10,10)
#'   SC=matrix(0,10,10)
#'   SC[locs]=temp[locs]
#'   diag(SC)=0
#' # Model fit
#'   fit<-SCFC(y,SC,method="QUIC",etaInd=1,siglam=10,sigmu=5,maxits=500,mu_init=0,a0_init=30,b0_init=6,c0=.5,outerits = 100);
#' @export

#Changes:
#transposing y
#rep(mu_init) -> mu
#etaind in conditions
#log(det()) -> determinant(logarithm=TRUE)



SCFC<- function(y,P,siglam=10,sigmu=5,maxits=500,method="glasso",etaInd=1,mu_init=NULL,a0_init=30,b0_init=5,c0=NULL,cov_init=NULL,alpha_init=NULL,outerits,eps=1e-4)
{
  y = t(y)
  Estimates=NULL
  T = nrow(y);
  p = ncol(y);
  niters=1;
  err=matrix(rep(NaN,maxits+1),maxits+1,1)
  err[1]=1;
  alpha=matrix(rep(0,p^2),ncol=p);
  endtime=matrix(rep(NaN,maxits+1),maxits+1,1);
  objfunc=matrix(rep(NaN,maxits+1),maxits+1,1);
  
  
  #Initialize covariance matrix,s, and inverse covariance matrix, OmegaN
  if(is.null(cov_init)){
    s=(1/T)*t(y)%*%y;#sample covariance matrix using observed functional data y?
  }else{
    s=cov_init
  }
  
  #Initialize mu
  if(is.null(mu_init)){
    mu=rep(0,p*(p-1)/2);
  }else{
    if(length(mu_init)==(p*(p-1)/2)){
      mu=mu_init
    } else{
      mu=rep(mu_init,p*(p-1)/2)
    }
  }
  
  #Initialize Omega
  hold_omeg=list(); hold_invomg=list();
  # grideta=seq(from=0.01,to=.8,length.out=25)
  grideta = seq(from=abs(min(s[upper.tri(s)]))/100,to=2*max(s[upper.tri(s)]),length.out = 25)
  fitBIC=rep(0,length(grideta))
  
  for(i in 1:length(grideta)){
    fitted=glasso(s,rho=grideta[i],penalize.diagonal=TRUE)
    hold_omeg[[i]]=fitted$wi
    hold_invomg[[i]]=fitted$w
    fitBIC[i]=-determinant(fitted$wi,logarithm = TRUE)$modulus+sum(diag(s%*%fitted$wi))+(log(T)/T)*sum(fitted$wi[upper.tri(fitted$wi)]!=0)
  }
  indloc=which(fitBIC==min(fitBIC))
  OmegaN=as.matrix(hold_omeg[[indloc]])
  invomg=as.matrix(hold_invomg[[indloc]])
  
  #Initialize alpha:
  alpha=log(grideta[indloc])*matrix(1,ncol=dim(P)[1],nrow=dim(P)[1])
  diag(alpha)=rep(-Inf,p)
  
  #Initialize eta
  if(etaInd){
    eta=stats::runif(1,min = .1,max = 5)
  }else {
    eta=0
  }
  
  
  #Initialize a0, b0
  if(is.null(a0_init)){
    a0=2
  }else {
    a0=a0_init}
  if(is.null(b0_init)){
    b0=2
  } else {
    b0=b0_init
  }
  if(is.null(c0)){
    cc=1
  }else{
    cc=c0
  }
  
  
  muN=mu;
  Psq=P*P;
  
  gam=1/(cc*mean(diag(OmegaN)))
  objfunc[niters]=-(T/2)*determinant(OmegaN,logarithm = TRUE)$modulus+(1/2)*sum(diag(OmegaN%*%s))+
    cc*exp(alpha)[upper.tri(exp(alpha))]%*%abs(OmegaN)[upper.tri(abs(OmegaN))] +
    (1/(2*siglam))*sum((alpha[upper.tri(alpha)]-muN+eta*P[upper.tri(P)])^2)
  -(a0-1)*log(eta)+b0*eta +(1/(2*sigmu))*sum((muN -mu)^2)-p*log(gam)+cc*gam*sum(diag(OmegaN))
  repeat {
    starttime=proc.time()[3]
    niters = niters + 1;
    
    ######
    Omegaold=OmegaN;
    alphaold=alpha;
    nuold=eta;
    muold=muN;
    
    
    
    ################################
    #Update--eta parameters
    ################################
    if(etaInd){
      AP=alphaold*P; #use updated alpha
      a=sum(Psq[upper.tri(Psq)])/(siglam)
      b=b0+(1/siglam)*(sum(AP[upper.tri(AP)])-muold%*%P[upper.tri(P)]);
      c=-(a0-1);
      eta=(-b+sqrt(b^2-4*a*c))/(2*a)
    }else {
      eta=0
    }
    ################################
    #Update--mu parameters
    ################################
    muN=( sigmu*(alphaold[upper.tri(alphaold)] + eta*P[upper.tri(P)]) + siglam*mu)/(siglam+sigmu)
    
    
    ##################################
    #Update--alpha parameters
    #    Goal:  update the upperdiagonal alpha's, then reform into pxp matrix (for GLasso procedure)
    ##################################
    timealp=proc.time()[3]
    alphaold=alpha;
    V=0*matrix(rep(1,p*p),ncol=p);
    omeg=Omegaold[upper.tri(Omegaold)];#use old Omega
    alpit=0;
    int=alphaold[upper.tri(alphaold)];
    
    while(alpit<1){# ofupdates to estimate of alpha before proceedig with next Omega estimation
      g=cc*siglam*abs(omeg)*exp(int) +(int-(muN -  eta*P[upper.tri(P)])) #first derivative
      H=1+cc*siglam*exp(int)*abs(omeg) #Hessian matrix:  diagonal matrix so only retain those values for computation gains
      dir=(1/H)*g
      
      ss=1
      alp_ele=(int-muN+eta*P[upper.tri(P)])
      f= 2*cc*siglam*sum(exp(int)*abs(omeg))+sum(alp_ele^2);
      repeat{
        #print(ss)
        nalpha=int-ss*dir;
        alp_ele_N=(nalpha-muN+eta*P[upper.tri(P)]);
        nf= 2*cc*siglam*sum(exp(nalpha)*abs(omeg))+sum((alp_ele_N)^2);
        if(f-nf>.45*ss*sum(g*dir) | sum(g^2)<1e-4 ){
          break
        }else{
          ss=.5*ss
          print(ss)
        }
        if(ss<.5^20){
          stop("too many ss iterations")
        }
      }
      alphaold=nalpha
      int=as.numeric(nalpha)
      alpit=alpit+1;
    }#end alpha while loop
    
    alpha1 = nalpha
    V[upper.tri(V)]=alpha1;
    V=V+t(V);
    alpha=V
    diag(alpha)=rep(log(2*(1/(gam))),p);
    
    
    
    
    ##################################
    #Update--inverse covariance matrix
    ##################################
    if(method=="QUIC"){
      OmegaNew=QUIC::QUIC(s,rho=(cc/2)*exp(alpha),X.init=Omegaold,W.init=invomg,maxIter=outerits,msg=0)
      OmegaN=OmegaNew$X#-->the estimated inverse covariance matrix
      invomg=OmegaNew$W
    }else if(method=="glasso"){
      OmegaNew=glasso::glasso(s ,rho = (cc/2)*exp(alpha),start ="warm",w.init=invomg,wi.init=Omegaold,penalize.diagonal = TRUE)
      OmegaN=OmegaNew$wi;
      invomg=OmegaNew$w;
    }
    gam= 1/(cc*mean(diag(OmegaN)));
    
    
    objfunc[niters]=-(T/2)*determinant(OmegaN,logarithm = TRUE)$modulus+.5*sum(diag(OmegaN%*%s))+
      cc*exp(alpha)[upper.tri(exp(alpha))]%*%abs(OmegaN)[upper.tri(abs(OmegaN))] +
      (1/(2*siglam))*sum((alpha[upper.tri(alpha)]-muN+eta*P[upper.tri(P)])^2)
    -(a0-1)*log(eta)+b0*eta +(1/(2*sigmu))*sum((muN -mu)^2)-p*log(gam)+cc*gam*sum(diag(OmegaN))
    
    endtime[niters-1]=proc.time()[3]-starttime;
    
    err[niters]=(objfunc[niters-1]-objfunc[niters])/abs(objfunc[niters])
    if(abs(err[niters])<eps || niters>maxits){
      
      if(abs(err[niters])<eps){
        Estimates$Method=method
        Estimates$iters=niters-1;
        Estimates$Omega=OmegaN;#precision matrix
        Estimates$Covariance=solve(OmegaN);#covariance matrix
        Estimates$lambda=exp(alpha);#shrinkage parameters
        Estimates$Eta=eta; #eta estimate:  average effect of struct on func
        Estimates$Mu=muN; #mu estimate:  controls overall sparcity
        Estimates$LogLike=objfunc[niters]; #Log-likelihood
        Estimates$del=abs(objfunc[niters-1]-objfunc[niters])/abs(objfunc[niters]); #change in obj func at covergence
        Estimates$Time=endtime[!is.na(endtime)]
        Estimates$Flag=0;
        Estimates$BIC=-determinant(OmegaN,logarithm = TRUE)$modulus+sum(s*OmegaN)+(log(T)/T)*sum(OmegaN[upper.tri(OmegaN)]!=0)
        break
      }
      else{
        message("Maximum Iterations exceeded")
        Estimates$Flag=1;
        break
      }
      #break
    }
    
    
  }#ends repeat loop
  
  return(Estimates)
}