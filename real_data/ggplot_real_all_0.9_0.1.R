library(ggplot2)

load("cqr-real-100_all-drp-wcp_0.9_0.1.rda")
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
dat2_50= do.call(data.frame,lapply(dat2, function(x) replace(x, is.infinite(x),50)))

load("real-split-cqr-100.rda")
len.mat=data[[2]]
up.all=data[[3]]
lo.all=data[[4]]

cov3 = data[[1]][,1]
len3 = data.frame(up.all[1,,]- lo.all[1,,])

full_drp_len = c(as.matrix(dat))
wcp_len = c(as.matrix(dat2_50))
split_drp_len = c(as.matrix(len3))

df.cov=data.frame( rep(cov.mat[,1],times = dim(up.all)[3]), rep(cov.mat[,2] ,times=dim(up.all)[3]), rep(cov3 ,times=dim(up.all)[3]))
df.len=data.frame( full_drp_len, wcp_len, split_drp_len)

bgnd <- theme_get()$panel.background$fill
names=c("Full DRP","WCP", "Split DRP")
colors_manual=c("red","dodgerblue4", "black")

#shape_manual=c(0,1,2,3,4,5)


col_names=c("cov.if","cov.wcp","cov.if3")
names(df.cov)=col_names
col_names=c("len.if","len.wcp","len.if3")
names(df.len)=col_names


nn1=nrow(df.len)
#df_names1 = data.frame(c(rep("synthetic data",R*2)), c(rep("Double robust prediction",R),rep("Weighted conformal prediction",R)))
df_names1 = data.frame(c(rep("Synthetic data",nn1)), c(rep(names[1],nn1)))
df_names2 = data.frame(c(rep("Synthetic data",nn1)), c(rep(names[2],nn1)))
df_names3 = data.frame(c(rep("Synthetic data",nn1)), c(rep(names[3],nn1)))



colnames(df_names1) = colnames(df_names2) = colnames(df_names3) = c("Setting", "Method")
df_names = rbind(df_names1,df_names2,df_names3)

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
rug_data$"V"=c(apply(df.cov,2,mean), apply(df.len,2,mean))
rug_data_cqr = rug_data

data_full=rbind(cov,len)
data_cqr = data_full

load("real-500-absolute_0.9_0.1.rda")
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
dat2_50= do.call(data.frame,lapply(dat2, function(x) replace(x, is.infinite(x),50)))

aa=rep(0,dim(dat)[2])
for (i in 1:dim(dat)[2]){
  aa[i]=sum(is.infinite(dat2[,i]))/dim(dat2)[1]
}
plot(aa)
mean(aa)

load("real-500-all-3splits.rda")
up.all=data[[3]]
lo.all=data[[4]]
cov3 = data[[1]][,1]
len3 = data.frame(up.all[1,,]- lo.all[1,,])

drp = wcp = wcp3 = c()
for (i in 1:length(dat)){
  drp=c(drp,dat[[i]])
  wcp=c(wcp,dat2_50[[i]])
  wcp3 = c(wcp3,len3[[i]])
}



df.cov2=data.frame( rep(cov.mat[,1],times = dim(up.all)[3]), rep(cov.mat[,2] ,times=dim(up.all)[3]), rep(cov3 ,times=dim(up.all)[3]) )
df.len2=data.frame( drp, wcp, wcp3)

df.cov = df.cov2
df.len = df.len2

bgnd <- theme_get()$panel.background$fill
names=c("Full DRP","WCP", "Split DRP")
colors_manual=c("red","dodgerblue4", "black")

#shape_manual=c(0,1,2,3,4,5)


col_names=c("cov.if","cov.wcp","cov.if3")
names(df.cov)=col_names
col_names=c("len.if","len.wcp","len.if3")
names(df.len)=col_names


nn2=nrow(df.len2)
#df_names1 = data.frame(c(rep("synthetic data",R*2)), c(rep("Double robust prediction",R),rep("Weighted conformal prediction",R)))
df_names1 = data.frame(c(rep("synthetic data",nn2)), c(rep(names[1],nn2)))
df_names2 = data.frame(c(rep("synthetic data",nn2)), c(rep(names[2],nn2)))
df_names3 = data.frame(c(rep("synthetic data",nn2)), c(rep(names[3],nn2)))



colnames(df_names1) = colnames(df_names2) = colnames(df_names3) = c("Setting", "Method")
df_names = rbind(df_names1,df_names2,df_names3)

cov_vec = as.vector(as.matrix(df.cov))
cov = cbind(cov_vec, df_names) 

len_vec = as.vector(as.matrix(df.len))
len = cbind(len_vec, df_names) 

cov$Var="Cov"
len$Var = "Width"

colnames(cov) = c("V","Setting","Method","Var")
colnames(len) = c("V","Setting","Method","Var")


var_names=c("Cov","Width")
settings = c("synthetic data")
rug_data = expand.grid(names,var_names,settings)
colnames(rug_data)=c("Method","Var","Setting")
rug_data$"V"=c(apply(df.cov2,2,mean), apply(df.len2,2,mean))


data_full=rbind(cov,len)
dummy_Cov <- data.frame(Var=c("Cov"),Z= 0.9)
data_absolute = data_full


rug_data_absolute = rug_data

data_full1 = data_absolute
data_full1$score = "Abs Residual"


n_factor = nrow(data_absolute)/nrow(data_cqr)
data_full2 = do.call("rbind", replicate(n_factor, data_cqr, simplify = FALSE))
data_full2$score = "CQR"
data_full = rbind(data_full1,data_full2)

rug_data_absolute$score = "Abs Residual"
rug_data_cqr$score = "CQR"
rug_full = rbind(rug_data_absolute, rug_data_cqr)

colors_manual2 = colors_manual[1:2]
final = ggplot(data_full, aes(x=V,fill = score)) +
  geom_histogram(aes(y=..count../sum(..count..)),position = "identity",alpha=0.5)+
  geom_rug(data=rug_data_absolute, sides="b",size=1,aes(x=V),color="red")+
  geom_rug(data=rug_data_cqr, sides="b",size=1,aes(x=V),color="dodgerblue4")+
  facet_grid(Method~Var,scales = "free",switch='x')+
  #scale_color_manual(values=colors_manual) +
  scale_fill_manual(values=colors_manual2)+
  ylab("Frequency")+
  theme(strip.text = element_text(size = 25, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        axis.title=element_text(size=25,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=20),
        axis.title.y = element_text(face="bold",size=25),
        axis.text.y = element_text(face="bold",size=20),
        #strip.text.x = element_blank(),
        legend.position="top",
        legend.title = element_text(size=30,face="bold"),
        legend.text = element_text(size=25),
        strip.background = element_blank(),
        strip.placement = "outside")+ 
  geom_vline(data=dummy_Cov, aes(xintercept=Z), linetype="dashed")
ggsave(("ggplot_real_all_0.9_0.1.pdf"),final, width = 16, height = 12, units = "in")
