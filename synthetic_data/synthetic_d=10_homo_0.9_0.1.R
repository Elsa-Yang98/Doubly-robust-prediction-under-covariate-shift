alpha = 0.05

method_r = c("SL.ridge")
#method_r = c("SL.bartMachine","SL.ridge","SL.glm")
method_pi = c("SL.mean","SL.randomForest","SL.glm")
method_m = c("SL.mean","SL.randomForest","SL.glm")
prop_train=0.3

Xfun <- function(n, d){
  matrix(runif(n * d), nrow = n, ncol = d)
}
psfun <- function(X){
  (1 + pbeta(X[, 1], 2, 4)) / 4
}
errdist <- rnorm
taufun <- function(X){
  2 / (1 + exp(-12 * (X[, 1] - 0.5))) * 2 / (1 + exp(-12 * (X[, 2] - 0.5)))
}
wsample = function(wts, ntest = length(wts)) {
  n = length(wts)
  i = c()
  while(length(i) < ntest) {
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

N = 1000
ntest = 10000
d = 10
dat.x=Xfun(2*N,d)
dat.y = taufun(dat.x) + errdist(N)

i = sample(2*N,N)
x = dat.x[i,]; y = dat.y[i]
x0 = dat.x[-i,]; y0 = dat.y[-i]

R = 100
cov.mat = len.mat = matrix(0,R,2)
up.all = lo.all = array(0,dim=c(2,R,ntest))

for (j in 1:R) {
  set.seed(j)
  cat("This is the", j,"out of", R,"runs")
  dat.x=Xfun(2*N,d)
  dat.y = taufun(dat.x) + errdist(N)
  
  i = sample(2*N,N)
  x = dat.x[i,]; y = dat.y[i]
  x0 = dat.x[-i,]; y0 = dat.y[-i]
  
  # Tilting
  i0 = wsample(psfun(x0)/(1-psfun(x0)), ntest = 10000)
  
  n0 = floor(length(i0))
  x00 = x0[i0,]; y00 = y0[i0]
  
  x_tmp = x_all = rbind(x,x00)
  t_all = as.factor(c(rep(0,nrow(x)),rep(1,nrow(x00))))
  
  ind_train_0 = sample(1:nrow(x),floor(nrow(x)*prop_train*2))
  ind_train_D1 = ind_train_0[1:floor(length(ind_train_0)/2)]
  ind_train_D20 = ind_train_0[-(1:floor(length(ind_train_0)/2))]
  ind_train_D21 =sample(1:nrow(x00),floor(nrow(x00)*prop_train))
  ind_train_D2 =c(ind_train_D20,nrow(x)+ind_train_D21)    
  
  
  x_train0=x[ind_train_0,]
  y_train0=y[ind_train_0]
  
  
  x_D20=x[ind_train_D20,]
  x_D2=x_all[ind_train_D2,]
  y_D2=y[ind_train_D20]
  t_D2=t_all[ind_train_D2]
  
  x_val = x_all[-c(ind_train_0,ind_train_D2),]
  
  x_val0 = x[-c(ind_train_0,ind_train_D20),]
  y_val0 = y[-c(ind_train_0,ind_train_D20)]
  t_val = t_all[-c(ind_train_0,ind_train_D2)]
  
  
  
  sl_r = SuperLearner(Y = as.numeric(y_train0), X = data.frame(x_train0),
                      SL.library = method_r)
  
  pred = predict(sl_r, data.frame(x_val0), onlySL = TRUE)
  
  r_D2 = abs(y_D2-predict(sl_r, data.frame(x_D20), onlySL = TRUE)$pred[,1])
  
  r_val = abs(y_val0-pred$pred[,1])
  
  r_0=c(r_D2,r_val)
  
  r_append_D3=c(r_val,rep(0,nrow(x00) - length(ind_train_D21) ))
  ########finish r#####
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
    prob = pmax(pmin(pred$pred[, 1], 0.9), 0.1)
    return(prob/(1-prob))
    #haha = w(xnew)
    # haha = rep(60,nrow(xnew))
    # return(haha)
  }
  
  
  
  
  
  f= function(theta){ mean(I(t_val==0)*pi_sl(x_val,as.numeric(t_D2)-1,x_D2,method=method_pi)*(I(r_append_D3<=theta)-m_sl(theta,x_val,as.numeric(r_D2<=theta),x_D20,method=method_m))+I(t_val==1)*(m_sl(theta,x_val,as.numeric(r_D2<=theta),x_D20,method=method_m)-(1-alpha)) )}
  #f= function(theta){ mean(I(t_val==0)*pi_np(x_val)*(I(r_val<=theta)-m(rep(theta, nrow(x_val)),x_val))+I(t_val==1)*(m(rep(theta, nrow(x_val)),x_val)-(1-alpha)) )}
  
  rhat = try(uniroot(f,c(1,max(r_0)/2),extendInt = "yes",tol=0.1)$root)
  if(class(rhat)!= "try-error"){
    rhat=rhat
  }else if (f(max(r_0)/2)>0){
    rhat = try(uniroot(f,c(0,4),tol=0.01)$root)
  }else{
    rhat = try(uniroot(f,c(max(r_0)/2,max(r_0)+0.1),tol=0.01)$root)
  } 
  
  
  
  pred_r_new = predict(sl_r, data.frame(x00), onlySL = TRUE)
  
  out_if_up=rhat+pred_r_new$pred[,1]
  out_if_lo=-rhat+pred_r_new$pred[,1]
  #out_if_up=rhat+predict(r_model,data.frame(x00))
  #out_if_lo=-rhat+predict(r_model,data.frame(x00))
  
  cov.mat[j,1] = mean(out_if_lo <= y00 & y00 <= out_if_up)
  len.mat[j,1] = median(out_if_up - out_if_lo)
  print(paste("j=",j,"DRP method, rhat=",rhat,"cov=",cov.mat[j,1],"len=",len.mat[j,1]))
  
  up.all[1,j,]=out_if_up
  lo.all[1,j,]=out_if_lo
  
  wts.glm = pi_sl(x_all,as.numeric(t_all)-1,x_all)
  #beta.mat1[j,] = coef(obj.glm)
  
  out = conformal.pred.sl.fast(x,x00, ind_train_0,m=1,pred=predict(sl_r, data.frame(x00), onlySL = TRUE)$pred[,1],res=r_val,alpha = 0.05, w=wts.glm)
  cov.mat[j,2] = mean(out$lo <= y00 & y00 <= out$up)
  len.mat[j,2] = median(out$up - out$lo)
  
  up.all[2,j,]=out$up
  lo.all[2,j,]=out$lo 
  
  print(paste("WCP method: cov=",cov.mat[j,2],"len=",len.mat[j,2]))
  
  
}

data=list(cov.mat[1:R,],len.mat[1:R,],up.all,lo.all)
save(data,file="synthetic_d=10_homo_0.9_0.1.rda")