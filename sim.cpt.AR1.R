sim.cpt.AR1=function(cpts,X,init,beta,sig2,nsim){
  # works for regression models with AR1 components
  
  # cpts must contain 0 and n
  # X is the design matrix for the regression part
  # init is the initial value for y from which the AR(1) is then calculated
  # beta is the regression fit (including the AR component as the last entry)
  # sig2 is the variance of the errors
  # nsim is the number of simulations
  
  cpts=cpts+1 # to take into account the initialization so same code can be used
  nseg=length(cpts)-1
  if(nseg==1){
    beta=matrix(beta,nrow=1)
  }
  if(is.null(dim(X))==T){
    p=1
  }
  else{
    p=ncol(X)    
  }
  if(length(sig2)==1){sig2=rep(sig2,nseg)}
  X=matrix(X,ncol=p) # makes sure it is a matrix so for a single regressor a vector can be given
  y=matrix(NA,ncol=nsim,nrow=cpts[length(cpts)]-1) #
  y[1,]=init # initialization same as original data
  seg=1
  for(i in 2:(cpts[length(cpts)]-1)){
    y[i,]=matrix(X[i-1,]%*%matrix(beta[seg,1:p],nrow=p),ncol=nsim,nrow=1)+beta[seg,p+1]*y[i-1,]+ rnorm(nsim,mean=0,sd=sqrt(sig2[seg]))
    if(i==cpts[seg+1]){seg=seg+1}
  }
  if(nsim==1){
    return(as.vector(y))
  }
  else{
    return(t(y))
  }
}
