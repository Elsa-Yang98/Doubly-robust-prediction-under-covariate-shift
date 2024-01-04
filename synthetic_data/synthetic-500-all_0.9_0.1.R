#source('conformal.pred.split.sl.R')
#estimate residuals in DRP and WCP efficiently using SuperLearner
#Use D1 to estimate r, whole data D1 and D2 to estimate two nuisance parameters pi and m, use D2 to evaluation f function
alpha=0.1
N=5000
n=3000
p=4

method_r = c("SL.ridge","SL.randomForest")
#method_r = c("SL.bartMachine","SL.ridge","SL.glm")
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



## Run the simulation
library(randomForest)
library(SuperLearner)

#To source, this is the modified function for the WCP method to train the residuals efficiently using SuperLearner
#In their original code they only use least squares to estimate the residuals
#utilize the prediction and the residuals fitted before
conformal.pred.sl.fast <- function (x,x0,ind_train0,m=1,pred,res, alpha=0.1,w=NULL){
  x = as.matrix(x)
  n = nrow(x)
  n0 = nrow(x0)
  i1=ind_train0
  i2 = (1:n)[-i1]
  res = matrix(res, ncol = m)
  pred = matrix(pred, ncol=m)
  # res = abs(y[i2] - matrix(predict.fun(out, x[i2, , drop = F]), nrow = n2))
  lo = up = matrix(0, n0, m)
  
  mad.x0 = rep(1, n0)
  
  weighted.quantile = function(v, prob, w=NULL, sorted=FALSE) {
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
  return(list(lo = lo, up = up))
}

conformal.pred.sl <- function (x, y, x0, ind_train0, alpha = 0.1, w = NULL,method=c("SL.bartMachine","SL.ridge","SL.glm")) 
{
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0, ncol = p)
  colnames(x0)=colnames(x)
  n0 = nrow(x0)
  if (is.null(w)) 
    w = rep(1, n + n0)
  
  
  i1 = ind_train0
  
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  
  
  sl_r =   SuperLearner(Y = as.numeric( y[i1] ), X = data.frame(x[i1, , drop = F]),
                        SL.library = method)
  fit = matrix(predict(sl_r, data.frame(x), onlySL = TRUE)$pred[,1],nrow=n)
  pred = matrix(predict(sl_r, data.frame(x0), onlySL = TRUE)$pred[,1],nrow=n0)
  
  res = abs(y[i2] - matrix(predict(sl_r, data.frame(x[i2, , drop = F]), onlySL = TRUE)$pred[,1], nrow = n2))
  
  
  # out = train.fun(x[i1, , drop = F], y[i1])
  # fit = matrix(predict.fun(out, x), nrow = n)
  # pred = matrix(predict.fun(out, x0), nrow = n0)
  m = ncol(pred)
  
  # res = abs(y[i2] - matrix(predict.fun(out, x[i2, , drop = F]), nrow = n2))
  lo = up = matrix(0, n0, m)
  
  mad.x0 = rep(1, n0)
  
  weighted.quantile = function(v, prob, w=NULL, sorted=FALSE) {
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
  return(list(pred = pred, lo = lo, up = up, fit = fit, split = i1))
}



conformal.pred.sl <- function (x, y, x0, ind_train0, alpha = 0.1, w = NULL,method=c("SL.bartMachine","SL.ridge","SL.glm")) 
{
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0, ncol = p)
  colnames(x0)=colnames(x)
  n0 = nrow(x0)
  if (is.null(w)) 
    w = rep(1, n + n0)
  
  
  i1 = ind_train0
  
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  
  
  sl_r =   SuperLearner(Y = as.numeric( y[i1] ), X = data.frame(x[i1, , drop = F]),
                        SL.library = method)
  fit = matrix(predict(sl_r, data.frame(x), onlySL = TRUE)$pred[,1],nrow=n)
  pred = matrix(predict(sl_r, data.frame(x0), onlySL = TRUE)$pred[,1],nrow=n0)
  
  res = abs(y[i2] - matrix(predict(sl_r, data.frame(x[i2, , drop = F]), onlySL = TRUE)$pred[,1], nrow = n2))
  
  
  # out = train.fun(x[i1, , drop = F], y[i1])
  # fit = matrix(predict.fun(out, x), nrow = n)
  # pred = matrix(predict.fun(out, x0), nrow = n0)
  m = ncol(pred)
  
  # res = abs(y[i2] - matrix(predict.fun(out, x[i2, , drop = F]), nrow = n2))
  lo = up = matrix(0, n0, m)
  
  mad.x0 = rep(1, n0)
  
  weighted.quantile = function(v, prob, w=NULL, sorted=FALSE) {
    if (is.null(w)) w = rep(1,length(v))
    if (!sorted) { o = order(v); v = v[o]; w = w[o] }
    i = which(cumsum(w/sum(w)) >= prob)
    if (length(i)==0) return(Inf) # Can happen with infinite weights
    else return(v[min(i)])
  }
  
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
  return(list(pred = pred, lo = lo, up = up, fit = fit, split = i1))
}





R = 500
cov.mat = len.mat = matrix(0,R,2)
beta.mat1 = beta.mat2 = matrix(0,R,p+1)
up.all = lo.all = array(0,dim=c(2,R,N-n))

beta = c(27.4,13.7,13.7,13.7)
for (j in 1:R) {
  set.seed(j)
  if (j %% 2 == 0) cat(j,".. ")
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
  
  # x_all[,1] = exp(x_tmp[,1] /2)
  # x_all[,2] = x_tmp[,2] /(1+exp(x_tmp[,2]))+10
  # x_all[,3] = (x_tmp[,1]*x_tmp[,3]/25+0.6)^3
  # x_all[,4] = (x_tmp[,2]+ x_tmp[,4] +20)^2
  # x = x_all[1:nrow(x),]
  # x00 = x_all[(1+nrow(x)):(nrow(x00)+nrow(x)),]
  
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
  
  sl_r = SuperLearner(Y = as.numeric(y_train0), X = data.frame(x_train0), SL.library = method_r)
  
  pred = predict(sl_r, data.frame(x_val0), onlySL = TRUE)
  
  r_train0 = abs(y_train0-predict(sl_r, data.frame(x_train0), onlySL = TRUE)$pred[,1])
  
  r_val = abs(y_val0-pred$pred[,1])
  r_0=rep(0,nrow(x))
  r_0[ind_train0]=r_train0
  r_0[-ind_train0]=r_val
  r_append=c(r_val,rep(0,nrow(x00) - length(ind_train1)))
  ####use linear regression to estimate r######
  
  
  ########finish r####
  
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
  

  
  
  pred_r_new = predict(sl_r, data.frame(x00), onlySL = TRUE)
  
 
  wts.glm = pi_sl(x_all,as.numeric(t_all)-1,x_all)
  #beta.mat1[j,] = coef(obj.glm)
  
  #out5 = conformal.pred.split.sl(x, y, x00, w=wts.glm,method=method_r)
  #out5 = conformal.pred.sl(x, y, x00, ind_train0,w=wts.glm,method=method_r)
  out = conformal.pred.sl.fast(x,x00, ind_train0,m=1,pred=predict(sl_r, data.frame(x00), onlySL = TRUE)$pred[,1],res=r_val,w=wts.glm)
  
  
  cov.mat[j,2] = mean(out$lo <= y00 & y00 <= out$up)
  len.mat[j,2] = median(out$up - out$lo)
  
  up.all[2,j,]=out$up
  lo.all[2,j,]=out$lo 
  
  print(paste("cov=",cov.mat[j,2],"len=",len.mat[j,2]))
  
  
}

data=list(cov.mat[1:R,],len.mat[1:R,],up.all,lo.all)
save(data,file="synthetic-500-all_0.9_0.1.rda")

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


