library(ggplot2)
library(ggpubr)
library(latex2exp)
load("CQR/data/synthetic-2000-split-cqr-conditional-100.rda")
set.seed(1234)
n=200
p=4
Xtest=matrix(rnorm(n*p),n,p)
xnorm = apply(Xtest,1,norm,type  ="2")
length(data)
dim(data[[1]])
cov_split = data[[1]][1,,,]
cov_split = apply(cov_split,3,mean)

cov_wcp = data[[1]][2,,,]
cov_wcp = apply(cov_wcp,3,mean)
n = length(cov_wcp)/2

n=200
cov_wcp = cov_wcp[1:n]
cov_split = cov_split[1:n]
cov_wcp_smooth = predict(smooth.spline(xnorm,cov_wcp),xnorm)$y
cov_split_smooth = predict(smooth.spline(xnorm,cov_split),xnorm)$y

len_split = data[[2]][1,,,] - data[[3]][1,,,]
len_wcp = data[[2]][2,,,] - data[[3]][2,,,]
len_split = apply(len_split,3,mean)
len_wcp = apply(len_wcp,3,mean)
len_split = len_split[1:n]
len_wcp_smooth = predict(smooth.spline(xnorm,len_wcp),xnorm)$y
len_split_smooth = predict(smooth.spline(xnorm,len_split),xnorm)$y
len_wcp = len_wcp[1:n]
len_mat= data.frame(length = c(len_split, len_wcp), ind = c(xnorm,xnorm), Method = c(rep("Split DRP", n), rep("WCP", n)))
len_mat_smooth = data.frame(length = c(len_split_smooth, len_wcp_smooth), ind = c(xnorm,xnorm), Method = c(rep("Split DRP", n), rep("WCP", n)))
p2_raw = ggplot(len_mat, aes(ind, length, linetype = Method, colour = Method))+
  geom_point()+
  geom_smooth(method = "loess")

p2 = ggplot(len_mat_smooth, aes(ind, length, linetype = Method, colour = Method))+
  geom_line(size = 1.2)
#geom_point(size = 5)
final2 = p2_raw+ scale_x_continuous(name=TeX(r'($l_2$ norm of test point)'), limits=c(0, 5), 
                                breaks = seq(0,5,by = 1),labels = seq(0,5,by = 1))+
  scale_y_continuous(name="Length", limits=c(30,110), 
                     breaks = seq(30,110,by = 20),labels = seq(30,110,by = 20))+
  theme_grey(base_size = 25)



xnorm = xnorm[1:n]

cov_mat = data.frame(coverage = c(cov_split, cov_wcp), ind = c(xnorm,xnorm), Method = c(rep("Split DRP", n), rep("WCP", n)))
cov_mat_smooth = data.frame(coverage = c(cov_split_smooth, cov_wcp_smooth), ind = c(xnorm,xnorm), Method = c(rep("Split DRP", n), rep("WCP", n)))
#ggsave(("ggplot_conditional_len.pdf"),final2, width = 20, height = 12, units = "in")
p_raw = ggplot(cov_mat, aes(x = ind, y = coverage, linetype = Method, colour = Method))+
  geom_point()+
  geom_smooth(method = "loess")+
  scale_x_continuous(name=TeX(r'($l_2$ norm of test point)'), limits=c(0, 5), 
                     breaks = seq(0,5,by = 1),labels = seq(0,5,by = 1))+
  theme(axis.text.x=element_text(colour=NA),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank()
  )

#geom_point(size = 3)
final_cov = p_raw+
  scale_y_continuous(name="Coverage", limits=c(0, 1), 
                     breaks = seq(0,1,by = 0.1),labels = seq(0,1,by = 0.1))+
  geom_hline(yintercept=0.9,linetype=2)+
  theme_grey(base_size = 25)+
  theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x = element_blank())

final_cond = ggarrange(final_cov, final2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
ggsave(("ggplot_conditional_smooth2.pdf"),final_cond, width = 20, height = 12, units = "in")


# yy2 = predict(smooth.spline(cov_split,xnorm),xnorm)$y
# plot(xnorm,yy2,ylim = c(0,1))
# yy3 = predict(smooth.spline(cov_split),cov_split)$y
# plot(xnorm,yy3,ylim = c(0.8,1))

cov_raw_smooth = ggplot(NULL, aes(x = ind, y = coverage, linetype = Method, colour = Method))+
  geom_point(data = cov_mat,size = 1.5)+
  geom_line(data = cov_mat_smooth, size = 1.5)+
  scale_x_continuous(name=TeX(r'($l_2$ norm of test point)'), limits=c(0, 5), 
                     breaks = seq(0,5,by = 1),labels = seq(0,5,by = 1))+
  scale_y_continuous(name="Coverage", limits=c(0, 1), 
                     breaks = seq(0,1,by = 0.1),labels = seq(0,1,by = 0.1))+
  geom_hline(yintercept=0.9,linetype=2)+
  theme(axis.text.x=element_text(colour=NA),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold",size=30),
        axis.title.y = element_text(face="bold",size=35),
        legend.position="top",
        legend.title = element_text(size=35,face="bold"),
        legend.text = element_text(size=30)
  )

len_raw_smooth = ggplot(NULL, aes(ind, length, linetype = Method, colour = Method))+
  geom_point(data = len_mat,size = 1.5)+
  geom_line(data = len_mat_smooth, size = 1.5)+
  scale_x_continuous(name=TeX(r'($L_2$ norm of test data X)'), limits=c(0, 5), 
                                    breaks = seq(0,5,by = 1),labels = seq(0,5,by = 1))+
  scale_y_continuous(name="Length", limits=c(30,110), 
                     breaks = seq(30,110,by = 20),labels = seq(30,110,by = 20))+
  theme(
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=30),
  axis.text.y = element_text(face="bold",size=30),
  legend.position="none",
  legend.title = element_text(size=35,face="bold"),
  legend.text = element_text(size=30)
)

final_cond_raw_smooth = ggarrange(cov_raw_smooth, len_raw_smooth, ncol=1, nrow=2)
ggsave(("ggplot_conditional_raw_smooth.pdf"),final_cond_raw_smooth, width = 16, height = 12, units = "in")


