#source('conformal.pred.split.sl.R')
#estimate residuals in DRP and WCP efficiently using SuperLearner
#Use D1 to estimate r, whole data D1 and D2 to estimate two nuisance parameters pi and m, use D2 to evaluation f function
rm(list = ls())
alpha=0.1
N=2000
n=600
p=4
trunc_list = c(0.9,1)
method_r = c("SL.ridge")
#method_r = c("SL.bartMachine","SL.ridge","SL.glm")
method_pi = c("SL.mean","SL.randomForest","SL.glm")
method_m = c("SL.mean","SL.randomForest","SL.glm")
prop_train=0.3

expit = function(u){
  exp(u)/(1+exp(u))
}

w = function(x,eta){
  #expit(x%*% c(-1,0.5,-0.25,-0.1))
  exp(x%*% c(-1,0.5,-0.25,-0.1)*eta)
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
    if(length(i)==1 || length(i)==0){
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





## Run the simulation
library(randomForest)
library(SuperLearner)
library(gbm)
#To source, this is the modified function for the WCP method to train the residuals efficiently using SuperLearner
#In their original code they only use least squares to estimate the residuals
#utilize the prediction and the residuals fitted before




R = 40
beta = c(27.4,13.7,13.7,13.7)

eta_list = seq(1,10,by = 1)
cov.mat = len.mat = array(0,dim=c(length(eta_list),length(trunc_list),R,2))
up.all = lo.all = cov.all = array(0,dim=c(length(eta_list),2,length(trunc_list),R,N-n))

truncation_mat = array(0,dim = c(length(eta_list), length(trunc_list),R))

for (l in 1:length(eta_list)){
  print(paste("this is the l =", l,'out of',length(eta_list), "th run"))
  eta = eta_list[l]
  for (j in 1:R) {
    set.seed(j)
    cat(j,".. ")
    dat.x=matrix(rnorm(N*p),N,p)
    dat.y = 210+dat.x%*%beta+rnorm(N)
    
    
    ind_split = sample(N,n)
    x = dat.x[ind_split,]; y = dat.y[ind_split]
    x0 = dat.x[-ind_split,]; y0 = dat.y[-ind_split]
    
    # Tilting
    i0 = wsample(w(x0,eta))
    
    
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
    
    
    
    pi_sl = function(xnew,ytrain,xtrain,method=method_pi, truncate_hi = truncate, truncate_lo = 0.1){
      #bw <- npcdistbw(r_train0~x_train0)
      xnew = data.frame(xnew)
      sl_pi = SuperLearner(Y = ytrain, X = data.frame(xtrain), family = binomial(),SL.library = method)
      pred = predict(sl_pi, xnew, onlySL = TRUE)
      prob = pred$pred[, 1]
      prob = pmax(pmin(pred$pred[, 1], truncate), truncate_lo)
      #prob = pmax(pmin(pred$pred[, 1], 1), )
      #return(prob)
      #prob = pmin(prob, truncate_hi)
      #prob = pmax(prob, truncate_lo)
      odds = (prob/(1-prob))
      return(pmin(odds, truncate_hi))
      #return(odds)
      #return(odds)
    }
    # df=data.frame(x_D2)
    # df$group = as.numeric(t_D2)-1
    # gbm1 <- gbm(group~., distribution ="bernoulli", data = df, cv.folds=5) 
    # pihat <- predict(gbm1, newdata = data.frame(x_D2))
    # pred = ifelse(pihat>0.5,1,0)
    # table(pred,df$group)
    # 
    # rf1 = randomForest(factor(as.numeric(t_D2)-1)~., data = data.frame(x_D2), proximity=T)
    # pihat2 = predict(rf1, data = data.frame(x_D2))
    # table(pihat2, df$group)
    # 
    t=0
    for (truncate in trunc_list){
      t=t+1
      tmp_pi = pi_sl(x_val,as.numeric(t_D2)-1,x_D2,method=method_pi, truncate_hi = truncate)
      #summary(pi_sl(x_val,as.numeric(t_D2)-1,x_D2,method=method_pi, truncate = 1))
      # sum(is.infinite(tmp))
      # summary((tmp))
      
      
      if (sum(is.infinite(tmp_pi))>0){
        rhat = max(r_append_D3)
        truncation_mat[l,t,j]=1
      }else{
        f= function(theta){ 
          tmp_m = m_sl(theta,x_val,as.numeric(r_D2<=theta),x_D20,method=method_m)
          mean(I(t_val==0)*tmp_pi*(I(r_append_D3<=theta)-tmp_m)+I(t_val==1)*(tmp_m-(1-alpha)) )}
        #f= function(theta){ mean(I(t_val==0)*pi_np(x_val)*(I(r_val<=theta)-m(rep(theta, nrow(x_val)),x_val))+I(t_val==1)*(m(rep(theta, nrow(x_val)),x_val)-(1-alpha)) )}
        
        rhat = try(uniroot(f,c(1,max(r_0)/2),extendInt = "yes",tol=0.1)$root)
        if(class(rhat)!= "try-error"){
          rhat=rhat
        }else if (f(max(r_0)/2)>0){
          tmp_f = try(f(0.01))
          if(class(tmp_f)!= "try-error"){rhat = try(uniroot(f,c(0.01,4),tol=0.01)$root)}
          else{rhat = try(uniroot(f,c(0.1,4),tol=0.01)$root)} 
        }else{
          rhat = try(uniroot(f,c(max(r_0)/2,max(r_0)+0.1),tol=0.01)$root)
        } 
      }
      
      
      
      
      pred_r_new = predict(sl_r, data.frame(x00), onlySL = TRUE)
      
      out_if_up=rhat+pred_r_new$pred[,1]
      out_if_lo=-rhat+pred_r_new$pred[,1]
      #out_if_up=rhat+predict(r_model,data.frame(x00))
      #out_if_lo=-rhat+predict(r_model,data.frame(x00))
      
      cov.mat[l,t,j,1] = mean(out_if_lo <= y00 & y00 <= out_if_up)
      len.mat[l,t,j,1] = median(out_if_up - out_if_lo)
      print(paste("DRP method: j=",j,"rhat=",rhat,"cov=",cov.mat[l,t,j,1],"len=",len.mat[l,t,j,1]))
      
      up.all[l,1,t,j,]=out_if_up
      lo.all[l,1,t,j,]=out_if_lo
      
      cov.all[l,1,t,j,] = out_if_lo <= y00 & y00 <= out_if_up
      
      
      #summary(pi_sl(x,as.numeric(t_all)-1,x_all,method=method_pi, truncate = truncate))
      #undebug(pi_sl)
      wts.glm = pi_sl(x_all,as.numeric(t_all)-1,x_all,method=method_pi, truncate_hi = truncate)
      #summary(wts.glm)
      #undebug(conformal.pred.sl.fast)
      #debug(conformal.pred.sl.fast)
      out = conformal.pred.sl.fast(x,x00, ind_train_0,m=1,pred=predict(sl_r, data.frame(x00), onlySL = TRUE)$pred[,1],res=r_val,alpha = alpha, w=wts.glm)
      #if (sum(is.infinite(out$up))>0){break}
      cov.mat[l,t,j,2] = mean(out$lo <= y00 & y00 <= out$up)
      len.mat[l,t,j,2] = median(out$up - out$lo)
      
      up.all[l,2,t,j,]=out$up
      lo.all[l,2,t,j,]=out$lo 
      
      cov.all[l,2,t,j,] = out$lo <= y00 & y00 <= out$up
      
      print(paste("WCP method: cov=",cov.mat[l,t,j, 2],"len=",len.mat[l,t,j, 2]))
    }
  }
}

data=list(cov.mat,len.mat,up.all,lo.all, cov.all, truncation_mat)
save(data,file="C:\\Users\\yacho\\Dropbox (Penn)\\truncation-propensity-upper-finite-0.9.rda")



