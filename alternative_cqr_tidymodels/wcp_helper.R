w = function(x) {
  # exp(x%*% c(-1,0.5,-0.25,-0.1))
  rep(1, nrow(x))
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

