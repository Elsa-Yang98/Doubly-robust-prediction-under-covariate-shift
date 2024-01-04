load("synthetic-500-all.rda")
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
dat2_50= do.call(data.frame,lapply(dat2, function(x) replace(x, is.infinite(x),10)))

load("synthetic-500-all-3splits.rda")
up.all=data[[3]]
lo.all=data[[4]]
cov3 = data[[1]][,1]
len3 = data.frame(up.all[1,,]- lo.all[1,,])

load("synthetic-500-all-2splits.rda")
up.all=data[[3]]
lo.all=data[[4]]
cov2 = data[[1]][,1]
len2 = data.frame(up.all[1,,]- lo.all[1,,])

drp = wcp = wcp3 = wcp2 = c()
for (i in 1:length(dat)){
  drp=c(drp,dat[[i]])
  wcp=c(wcp,dat2_50[[i]])
  wcp3 = c(wcp3,len3[[i]])
  wcp2 = c(wcp2,len2[[i]])
}



df.cov2=data.frame( rep(cov.mat[,1],times = dim(up.all)[3]), rep(cov.mat[,2] ,times=dim(up.all)[3]), rep(cov3 ,times=dim(up.all)[3]), rep(cov2 ,times=dim(up.all)[3]) )
df.len2=data.frame( drp, wcp, wcp3, wcp2)





df.cov=rbind(df.cov2)
df.len=rbind(df.len2)


bgnd <- theme_get()$panel.background$fill
names=c("DRP w. full data","WCP", "DRP w. three splits","DRP w. two splits" )


#shape_manual=c(0,1,2,3,4,5)


col_names=c("cov.if","cov.wcp","cov.if3")
names(df.cov)=names
col_names=c("len.if","len.wcp","len.if3")
names(df.len)=names

settings = c("synthetic data")
#df_names = expand.grid(settings,names)

nn2=nrow(df.len2)
#df_names1 = data.frame(c(rep("synthetic data",R*2)), c(rep("Double robust prediction",R),rep("Weighted conformal prediction",R)))
df_names1 = data.frame(rep("synthetic data",nn2), c(rep(names[1],nn2)))
df_names2 = data.frame(rep("synthetic data",nn2), c(rep(names[2],nn2)))
df_names3 = data.frame(rep("synthetic data",nn2), c(rep(names[3],nn2)))
df_names4 = data.frame(rep("synthetic data",nn2), c(rep(names[4],nn2)))

colnames(df_names1) = colnames(df_names2) = colnames(df_names3) =  colnames(df_names4) = c("Setting", "Method")
df_names = rbind(df_names1,df_names2,df_names3, df_names4)

cov_vec = as.vector(as.matrix(df.cov))
cov = cbind(cov_vec, df_names) 

len_vec = as.vector(as.matrix(df.len))
len = cbind(len_vec, df_names) 

cov$Var="Coverage"
len$Var = "Width"

colnames(cov) = c("V","Setting","Method","Var")
colnames(len) = c("V","Setting","Method","Var")


var_names=c("Coverage","Width")
rug_data = expand.grid(names,var_names,settings)
colnames(rug_data)=c("Method","Var","Setting")
rug_data$"V"=c(apply(df.cov2,2,mean), apply(df.len2,2,mean))


data_full=rbind(cov,len)
dummy_coverage <- data.frame(Var=c("Coverage"),Z= 0.9)

colors_manual=c("purple", "red","dodgerblue4","black")
final_plot_synthetic = ggplot(data=data_full, aes(x=V,group=Method,color=Method,fill=Method)) +
  geom_histogram(aes(y=..count../sum(..count..)),alpha=0.5)+
  geom_rug(data=rug_data, sides="b",size=1,aes(x=V,group=Method,color=Method))+
  facet_grid(Method~Var,scales = "free",switch='x')+
  scale_color_manual(values=colors_manual) +
  scale_fill_manual(values=colors_manual)+
  ylab("Frequency")+
  theme(strip.text = element_text(size = 18, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        axis.title=element_text(size=18,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=12),
        #strip.text.x = element_blank(),
        legend.position="top",
        legend.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside")+ 
  geom_vline(data=dummy_coverage, aes(xintercept=Z), linetype="dashed")

ggsave(("ggplot_synthetic2_new.pdf"),final_plot_synthetic, width = 16, height = 12, units = "in")


