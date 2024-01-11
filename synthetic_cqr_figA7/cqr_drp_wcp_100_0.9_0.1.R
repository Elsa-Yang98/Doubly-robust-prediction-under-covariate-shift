rm(list = ls())
library(randomForest)
library("conformalInference")
library(SuperLearner)
alpha=0.1 # target miscoverage level
N=5000
n=3000
p=4

method_pi = c("SL.mean","SL.randomForest","SL.glm")
method_m = c("SL.mean","SL.randomForest","SL.glm")
prop_train=0.5

expit = function(u){
  exp(u)/(1+exp(u))
}

w = function(x) {
  #expit(x%*% c(-1,0.5,-0.25,-0.1))
  exp(x%*% c(-1,0.5,-0.25,-0.1))
}


wsample = function(wts, frac=1) {
  n = length(wts)
  i = c()
  while(length(i) < n*frac) {
    j = which( runif(n) <= wts/sum(wts)  )[1] 
    if ( !is.na(j) ) {
      i = c(i, j )
    }
  }
  return(i)
}


#This is the modified function for the WCP method to train the residuals efficiently using SuperLearner
#In their original code they only use least squares to estimate the residuals
#utilize the prediction and the residuals fitted before.
conformal.pred.sl.fast <- function (x,x0,ind_train0,m=1,pred,res, alpha=0.1,w=NULL,type="mean", quantiles = quantiles){
  # based on tibshiraniri's code, type = either "CQR" or "mean", pred are the point predictions, res are the scores that we will be finding the weighted quantile of
  x = as.matrix(x)
  n = nrow(x)
  n0 = nrow(x0)
  i1=ind_train0
  i2 = (1:n)[-i1]
  res = matrix(res, ncol = m)
  
  # res = abs(y[i2] - matrix(predict.fun(out, x[i2, , drop = F]), nrow = n2))
  lo = up = matrix(0, n0, m)
  
  mad.x0 = rep(1, n0)
  
  weighted.quantile = function(v, prob, w=NULL, sorted=FALSE){
    if (is.null(w)) w = rep(1,length(v))
    if (!sorted) { o = order(v); v = v[o]; w = w[o] }
    t = pmax(pmin(w, 10), 0.05)
    ww = t/sum(t)
    i = which(cumsum(ww) >= prob)
    #if (length(i)==0) return(Inf) # Can happen with infinite weights
    if(length(i)==1){
      return(max(v[-length(v)]))
    }else{
      return(v[min(i)])
    } 
  }
  
  
  if (type == "mean"){
    pred = matrix(pred, ncol=m)
    for (l in 1:m) {
      
      o = order(res[, l])
      r = res[o, l]
      ww = w[i2][o]
      for (i in 1:n0) {
        q = weighted.quantile(c(r, Inf), 1 - alpha, w = c(ww, w[n + i]), sorted = TRUE)
        lo[i, l] = pred[i, l] - q * mad.x0[i]
        up[i, l] = pred[i, l] + q * mad.x0[i]
      }
    }
  }
  
  if (type == "CQR"){
    if (length(quantiles)!=2){print("error: number of quantiles should be 2")}
    if (m==1){
      l=1
      pred = matrix(pred,ncol = 2)
      o = order(res[, l])
      r = res[o, l]
      ww = w[i2][o]
      for (i in 1:n0) {
        q = weighted.quantile(c(r, Inf), 1 - alpha, w = c(ww, w[n + i]), sorted = TRUE)
        lo[i, l] = pred[i, 1] - q * mad.x0[i]
        up[i, l] = pred[i, 2] + q * mad.x0[i]
      }
    }
  }
  else{
    print("error: multidimensional cqr")
  }
  return(list(lo = lo, up = up))
}


R = 100
cov.mat = len.mat = matrix(0,R,2)
up.all = lo.all = array(0,dim=c(2,R,N-n))

beta = c(27.4,13.7,13.7,13.7)
for (j in 1:R) {
  set.seed(j)
  cat(j,".. ")
  dat.x=matrix(rnorm(N*p),N,p)
  dat.y = 210+dat.x%*%beta+rnorm(N)
  
  
  i = sample(N,n)
  x = dat.x[i,]; y = dat.y[i]
  x0 = dat.x[-i,]; y0 = dat.y[-i]
  
  # Tilting
  i0 = wsample(w(x0))
  
  #i0 = wsample(w(dat.z[-i,]))
  n0 = floor(length(i0))
  x00 = x0[i0,]; y00 = y0[i0]
  
  x_tmp = x_all = rbind(x,x00)

  
  t_all = as.factor(c(rep(0,nrow(x)),rep(1,nrow(x00))))
  
  
  ind_train0=sample(1:nrow(x),floor(nrow(x)*prop_train))
  ind_train1=sample(1:nrow(x00),floor(nrow(x00)*prop_train))
  ind_train =c(ind_train0,nrow(x)+ind_train1)    
  
  t_train=t_all[ind_train]
  x_train=x_all[ind_train,]
  
  x_train0=x[ind_train0,]
  y_train0=y[ind_train0]
  
  x_val = x_all[-ind_train,]
  
  x_val0 = x[-ind_train0,]
  y_val0 = y[-ind_train0]
  t_val = t_all[-ind_train]
  
  #r=abs(y - x[,1])   #assume the simplest model for r, an arbitrary function of x and y
  ####use linear regression to estimate r######
  
  y_train0=y[ind_train0]
  
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
  fit_r = grf::quantile_forest(x_train0,y_train0, quantiles = quantiles)
  
  
  quantile_train0 = predict(fit_r, x_train0, quantiles = quantiles)$predictions
  r_train0 = pmax( quantile_train0[,1]- y_train0, y_train0- quantile_train0[,2] )
  
  quantile_val = predict(fit_r, x_val0, quantiles = quantiles)$predictions
  r_val = pmax( quantile_val[,1]- y_val0, y_val0- quantile_val[,2] )
  
  r_0 = rep(0,nrow(x))
  r_0[ind_train0] = r_train0
  r_0[-ind_train0] = r_val
  r_append=c(r_val,rep(0,nrow(x00) - length(ind_train1)))
  ########finish the scores r####
  # estimate two nuisance functions 
  m_sl=function(theta,xnew,ytrain,xtrain,method=method_m){
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
    prob = pmax(pmin(prob, 0.9), 0.1)
    return(prob/(1-prob))
  }
  
  f= function(theta){ mean(I(t_val==0)*pi_sl(x_val,as.numeric(t_all)-1,x_all)*(I(r_append<=theta)-m_sl(theta,x_val,as.numeric(r_0<=theta),x))+I(t_val==1)*(m_sl(theta,x_val,as.numeric(r_0<=theta),x)-(1-alpha)) )}
  #f= function(theta){ mean(I(t_val==0)*pi_np(x_val)*(I(r_val<=theta)-m(rep(theta, nrow(x_val)),x_val))+I(t_val==1)*(m(rep(theta, nrow(x_val)),x_val)-(1-alpha)) )}
  
  rhat = try(uniroot(f,c(sort(r_0)[8],max(r_0)/2),extendInt = "yes",tol=0.1)$root)
  if(class(rhat)!= "try-error"){
    rhat=rhat
  }else if (f(max(r_0)/2)>0){
    rhat = try(uniroot(f,c(0,4),tol=0.01)$root)
  }else{
    rhat = try(uniroot(f,c(max(r_0)/2,max(r_0)+0.1),tol=0.01)$root)
  } 

  pred_r_new = predict(fit_r, data.frame(x00), quantiles = quantiles)$predictions
  
  out_if_up=rhat+pred_r_new[,2]
  out_if_lo=-rhat+pred_r_new[,1]

  cov.mat[j,1] = mean(out_if_lo <= y00 & y00 <= out_if_up)
  len.mat[j,1] = median(out_if_up - out_if_lo)
  print(paste("j=",j,"DRP produces rhat=",rhat,"cov=",cov.mat[j,1],"len=",len.mat[j,1]))
  
  up.all[1,j,]=out_if_up
  lo.all[1,j,]=out_if_lo
  
  wts.glm = pi_sl(x_all,as.numeric(t_all)-1,x_all)

  out = conformal.pred.sl.fast(x,x00, ind_train0,m=1,pred=pred_r_new, res=r_val,w=wts.glm, type ="CQR", quantiles = quantiles)
  
  cov.mat[j,2] = mean(out$lo <= y00 & y00 <= out$up)
  len.mat[j,2] = median(out$up - out$lo)
  
  up.all[2,j,]=out$up
  lo.all[2,j,]=out$lo 
  
  print(paste("WCP: cov=",cov.mat[j,2],"len=",len.mat[j,2]))
  
  
}

data=list(cov.mat[1:R,],len.mat[1:R,],up.all,lo.all)
save(data,file="cqr-synthetic-100-drp-wcp_0.9_0.1.rda")

RR=R
A2=len.mat[1:RR,1]
B2=len.mat[1:RR,2]
A=cov.mat[1:RR,1]
B=cov.mat[1:RR,2]

mean(A)
mean(B)
sd(A)
sd(B)
mean(A2)
mean(B2)
sd(A2)
sd(B2)


