rm(list = ls())
load('truncation-propensity-upper-infinite.rda')
library(tidyverse)
library(scales)
library(latex2exp)

cov.mat = data[[1]]
len.mat = data[[2]]
up.all = data[[3]]
lo.all = data[[4]]
cov.all = data[[5]]
truncation_mat = data[[6]]
eta_list = seq(10)

cov.mat = cov.mat[1:9,,,]
len.mat = len.mat[1:9,,,]
up.all = up.all[1:9,,,,]
lo.all = lo.all[1:9,,,,]
cov.all = cov.all[1:9,,,,]
eta_list = eta_list[1:9]


load('truncation-propensity-upper-finite-0.9.rda')

cov.mat2 = data[[1]]
len.mat2 = data[[2]]
up.all2 = data[[3]]
lo.all2 = data[[4]]
cov.all2 = data[[5]]
truncation_mat = data[[6]]
eta_list = seq(9)

cov.mat2 = cov.mat2[1:9,,,]
len.mat2 = len.mat2[1:9,,,]
up.all2 = up.all2[1:9,,,,]
lo.all2 = lo.all2[1:9,,,,]
cov.all2 = cov.all2[1:9,,,,]
eta_list2 = eta_list[1:9]

means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  b[!is.finite(b)] = 10
  rowMeans(b, dims = n - 1)
}
library(tibble)
cov_drp = data.frame("eta" = seq(length(eta_list)),means.along(cov.mat[,,,1], 3)[,c(1)], "Method" = "DRP", "Type" = "Coverage")
cov_drp = add_column(cov_drp, "0.9" = means.along(cov.mat2[,,,1], 3)[,c(1)], .after = 2)

cov_wcp = data.frame("eta" = seq(length(eta_list)),means.along(cov.mat[,,,2], 3)[,c(1)], "Method" = "WCP", "Type" = "Coverage")
cov_wcp = add_column(cov_wcp, "0.9" = means.along(cov.mat2[,,,2], 3)[,c(1)], .after = 2)

len_drp = data.frame("eta" = seq(length(eta_list)),means.along(len.mat[,,,1], 3)[,c(1)], "Method" = "DRP", "Type" = "Length")
len_drp = add_column(len_drp, "0.9" = means.along(len.mat2[,,,1], 3)[,c(1)], .after = 2)

len_wcp = data.frame("eta" = seq(length(eta_list)),means.along(len.mat[,,,2], 3)[,c(1)], "Method" = "WCP", "Type" = "Length")
len_wcp = add_column(len_wcp, "0.9" = means.along(len.mat2[,,,2], 3)[,c(1)], .after = 2)

names(cov_drp)[2:3] = c("0.95", "0.9")
names(cov_wcp)[2:3] = c("0.95", "0.9")
names(len_drp)[2:3] = c("0.95", "0.9")
names(len_wcp)[2:3] = c("0.95", "0.9")
tmp = rbind(cov_drp, cov_wcp, len_drp, len_wcp)
names(tmp)[2:3] = c("0.95", "0.9")

df1 = tmp %>%
  #select(alpha, forster_poly, ls_poly,forster_ns,ls_ns, forster_bs, ls_bs) %>%
  gather("Truncation","Value",2:3, -eta)


list.ab <- c()

# DOUBLE LOOP TO CREATE RECURSIVE COMBINATIONS USING append
ct=0    
for( i in 1:length(df1['Method'][[1]]) ) {    
  ct=ct+1
  list.ab[ct] <- paste(df1['Method'][[1]][[i]], df1['Truncation'][[1]][[i]])     
}
df1['type'] = list.ab
dummy_Cov <- data.frame(Var=c("Cov"),Z= 0.90)

final = ggplot(df1, aes(x = eta, y = Value)) +
  geom_line(aes(color = type, linetype = type), linewidth = 2)+
  scale_linetype_manual(values = c(rep("solid",2), rep("dashed", 2)))+
  #scale_color_manual(values=c(rep(hue_pal()(11)[1:11],2)))+
  facet_grid(Type~.,scales = "free",switch='y')+
  geom_point(aes(color = type, shape = type), size = 5)+
  #ylab("n*MSE")+
  xlab(TeX(r'( $\eta$ in the propensity score)'))+
  theme(strip.text = element_text(size = 40, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        legend.position="top",
        legend.title = element_text(size=30,face="bold"),
        legend.text = element_text(size=30),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x=element_text(size=35,face="bold", vjust = -0.6),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=35),
        axis.text.y = element_text(face="bold",size=35),
        #axis.title.y =  element_text(face="bold",size=40),
        axis.title.y =  element_blank(),
        #strip.text.x = element_blank(),
        legend.key.size = unit(5,"line")
  )+
  geom_hline(data=df1 %>% filter(Type == "Coverage"), aes(yintercept=0.9), linetype="dashed")+
  scale_x_continuous(breaks=seq(10))


ggsave(("ggplot_truncation_propensity_finite.pdf"),final, width = 16, height = 12, units = "in")


