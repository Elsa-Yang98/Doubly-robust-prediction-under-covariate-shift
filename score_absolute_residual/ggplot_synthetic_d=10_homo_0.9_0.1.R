rm(list = ls())
library(ggplot2)
load("synthetic_d=10_homo_quantile_only.rda")
cov_q = data[[1]]
len_q = data[[2]]
up_q = data[[3]]
lo_q = data[[4]]


load('synthetic_d=10_homo_0.9_0.1.rda')
cov.mat = data[[1]]
len.mat=data[[2]]
up.all=data[[3]]
lo.all=data[[4]]
R=dim(cov.mat)[1]

dat= up.all[1,,]- lo.all[1,,]
dat=data.frame(dat)
dat1= do.call(data.frame,lapply(dat, function(x) replace(x, is.infinite(x),10)))
dat2= up.all[2,,]- lo.all[2,,]
dat2=data.frame(dat2)
dat2_10= do.call(data.frame,lapply(dat2, function(x) replace(x, is.infinite(x),10)))

dat3 = data.frame(up_q - lo_q)

drp = wcp = quantile = c()
for (i in 1:length(dat)){
  drp=c(drp,dat[[i]])
  wcp=c(wcp,dat2_10[[i]])
}
quantile = unlist(dat3)

df.cov0=data.frame( rep(cov.mat[,1],times = dim(up.all)[3]), rep(cov.mat[,2] ,times=dim(up.all)[3]),  rep(cov_q ,times=dim(up.all)[3]))
df.len0=data.frame( drp, wcp, quantile)

df.cov = df.cov0
df.len = df.len0

bgnd <- theme_get()$panel.background$fill
names=c("Split DRP","WCP","Quantile only")
#colors_manual=c("red","dodgerblue4")

#shape_manual=c(0,1,2,3,4,5)


col_names=c("cov.if","cov.wcp", 'cov.q')
names(df.cov)=col_names
col_names=c("len.if","len.wcp","len.q")
names(df.len)=col_names


nn1=nrow(df.len0)
#df_names1 = data.frame(c(rep("synthetic data",R*2)), c(rep("Double robust prediction",R),rep("Weighted conformal prediction",R)))
df_names1 = data.frame(c(rep("Synthetic data",nn1)), c(rep(names[1],nn1)))
df_names2 = data.frame(c(rep("Synthetic data",nn1)), c(rep(names[2],nn1)))
df_names3 = data.frame(c(rep("Synthetic data",nn1)), c(rep(names[3],nn1)))


colnames(df_names1) = colnames(df_names2) = colnames(df_names3) = c("Setting", "Method")
df_names = rbind(df_names1,df_names2, df_names3)

cov_vec = as.vector(as.matrix(df.cov))
cov = cbind(cov_vec, df_names) 

len_vec = as.vector(as.matrix(df.len))
len = cbind(len_vec, df_names) 

cov$Var="Cov"
len$Var = "Width"

colnames(cov) = c("V","Setting","Method","Var")
colnames(len) = c("V","Setting","Method","Var")


var_names=c("Cov","Width")
settings = c("Synthetic data")
rug_data = expand.grid(names,var_names,settings)
colnames(rug_data)=c("Method","Var","Setting")
rug_data$"V"=c(apply(df.cov0,2,mean), apply(df.len0,2,mean))
rug_data_absolute = rug_data

data_full=rbind(cov,len)
data_absolute = data_full

data_full1 = data_absolute
data_full1$score = "Abs Residual"

data_full = data_full1




rug_data_absolute$score = "Abs Residual"
rug_full = rug_data_absolute
dummy_Cov <- data.frame(Var=c("Cov"),Z= 0.95)


final = ggplot(data=data_full, aes(x=V,group=Method,color=Method,fill=Method)) +
  geom_histogram(aes(y=..count../sum(..count..)),alpha=0.5)+
  geom_rug(data=rug_data_absolute, sides="b",size=1,aes(x=V,group=Method,color=Method))+
  facet_grid(Method~Var,scales = "free",switch='x')+
  #scale_color_manual(values=colors_manual) +
  #scale_fill_manual(values=colors_manual)+
  ylab("Frequency")+
  theme(strip.text = element_text(size = 22, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        axis.title=element_text(size=25,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        #strip.text.x = element_blank(),
        legend.position="top",
        legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside")+ 
  geom_vline(data=dummy_Cov, aes(xintercept=Z), linetype="dashed")
ggsave(("ggplot_synthetic_d=10_homo_3methods_0.9_0.1.pdf"),final, width = 16, height = 12, units = "in")

