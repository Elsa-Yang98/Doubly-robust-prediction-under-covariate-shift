library(cfcausal)
library(SuperLearner)

set.seed(1234)
alpha=0.1
N=2000
n=200
p=4
R = 100
beta = c(27.4,13.7,13.7,13.7)
cov.all = up.all = lo.all = array(0,dim=c(2,R,R,n))

method_pi = c("SL.mean","SL.randomForest","SL.glm")
method_m = c("SL.mean","SL.randomForest","SL.glm")


expit = function(u){
  exp(u)/(1+exp(u))
}
j=1

Xtest=matrix(rnorm(n*p),n,p)



for (i in 1:(2*R)){
  print(i)
  X=matrix(rnorm(0.4*N*p),0.4*N,p)
  Y = 210+X%*%beta+rnorm(0.4*N)
  
  ps = expit(X%*% c(-1,0.5,-0.25,-0.1)) #missing probability
  T <- as.numeric(runif(0.4*N)<ps)
  Y[T==1] <- NA
  
  quantRF <- function(Y, X, Xtest, quantiles, ...){
    fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
    res <- predict(fit, Xtest, quantiles = quantiles)$predictions
    if (length(quantiles) == 1){
      res <- as.numeric(res)
    } else {
      res <- as.matrix(res)
    }
    return(res)
  }
  
  quantiles = c(0.05, 0.95)
  fit_r = grf::quantile_forest(X[!is.na(Y),], Y[!is.na(Y)], quantiles = quantiles)
  res_r = predict(fit_r, Xtest, quantiles = quantiles)$predictions
  
  quantile_observe = predict(fit_r, X[!is.na(Y),], quantiles = quantiles)$predictions
  R_observe = pmax( quantile_observe[,1]- Y[!is.na(Y)], Y[!is.na(Y)] - quantile_observe[,2] )
  R_all = rep(0,0.4*N)
  R_all[!is.na(Y)] = R_observe
  
  
  X_train=matrix(rnorm(0.6*N*p),0.6*N,p)
  Y_train = 210+X_train%*%beta+rnorm(0.6*N)
  ps_train = expit(X_train%*% c(-1,0.5,-0.25,-0.1)) #missing probability
  T_train <- as.numeric(runif(0.6*N)<ps_train)
  Y_train[T_train==1] <- NA
  
  fit_r_train = grf::quantile_forest(X_train[!is.na(Y_train),], Y_train[!is.na(Y_train)], quantiles = quantiles)
  res_r_train = predict(fit_r_train, Xtest, quantiles = quantiles)$predictions
  
  quantile_observe_train = predict(fit_r_train, X_train[!is.na(Y_train),], quantiles = quantiles)$predictions
  R_observe_train = pmax( quantile_observe_train[,1]- Y_train[!is.na(Y_train)], Y_train[!is.na(Y_train)] - quantile_observe_train[,2] )
  R_all_train = rep(0,0.6*N)
  R_all_train[!is.na(Y_train)] = R_observe_train

  
  m_sl=function(xnew,ytrain,xtrain,method=method_m){
    xnew = data.frame(xnew)
    sl_m = SuperLearner(Y = ytrain, X = data.frame(xtrain), family = binomial(),
                        SL.library = method)
    pred = predict(sl_m, xnew, onlySL = TRUE)
    return(pred$pred[,1])
  }
  
  pi_sl = function(xnew,ytrain,xtrain,method=method_pi){
    #bw <- npcdistbw(r_train0~x_train0)
    xnew = data.frame(xnew)
    sl_pi = SuperLearner(Y = ytrain, X = data.frame(xtrain), family = binomial(),SL.library = method)
    pred = predict(sl_pi, xnew, onlySL = TRUE)
    prob = pred$pred[, 1]
    prob = pmax(pmin(prob, 0.99), 0.01)
    return(prob/(1-prob))
  }
  
  
  
  f= function(theta){ mean(I(T_train==0)*pi_sl(X_train,T,X)*(I(R_all_train <=theta)-m_sl(X_train,as.numeric(R_all<=theta),X))+I(T_train==1)*(m_sl(X_train,as.numeric(R_all<=theta),X)-(1-alpha)) )}
  
  rhat = try(uniroot(f,c(sort(R_all_train)[8],max(R_all_train)-0.1),extendInt = "yes",tol=0.1)$root)
  if(class(rhat)!= "try-error"){
    rhat=rhat
  }else if (f(max(R_all)-0.1)>0){
    rhat = try(uniroot(f,c(min(R_all_train)+0.1,sort(R_all_train)[8]),tol=0.01)$root)
  }else{
    rhat = try(uniroot(f,c(max(R_all_train)-0.1,max(R_all_train)),tol=0.01)$root)
  } 
  
  if(class(rhat)!= "try-error"){
    j=j+1
    if (j==R+1){break}
  }else{ next}
  

for (k in 1:R){
  Ytest = 210+Xtest%*%beta+rnorm(n)
  
  out_if_up=rhat+res_r[,2]
  out_if_lo=-rhat+res_r[,1]
  
  up.all[1,j,k,]=out_if_up
  lo.all[1,j,k,]=out_if_lo
  
  cov.all[1,j,k,] = out_if_lo  <= Ytest & Ytest <= out_if_up
  print(paste("j=",j,"rhat=",rhat))
  
  obj <- conformalCf(X, Y, estimand="missing", type = "CQR", 
                     quantiles = c(0.05, 0.95),
                     outfun = "quantRF", useCV = FALSE)
  wcp_res = predict(obj, Xtest, alpha = 0.1)
  up.all[2,j,k,]= wcp_res[,2]
  lo.all[2,j,k,]= wcp_res[,1]
  cov.all[2,j,k,] = wcp_res[,1] <= Ytest & Ytest <= wcp_res[,2]
}
}

data=list(cov.all,up.all,lo.all)
save(data,file="synthetic-2000-split-cqr-conditional-100.rda")