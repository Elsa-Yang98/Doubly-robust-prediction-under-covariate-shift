source('wcp_helper.R')
library(tidymodels)
library(probably)
# Modified based on https://www.tidymodels.org/learn/models/conformal-regression/
coverage <- function(x) {
  x %>% 
    mutate(in_bound = .pred_lower <= y & .pred_upper >= y) %>% 
    summarise(coverage = mean(in_bound) * 100)
}
width <- function(x) {
  x %>% 
    mutate(size = .pred_upper - .pred_lower) %>% 
    summarise(width = mean(size))
}
N=6000
n=3000
p=4
beta = c(27.4,13.7,13.7,13.7)
dat.x=matrix(rnorm(N*p),N,p)
dat.y = 210+dat.x%*%beta+rnorm(N)
prop_train = 0.5

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

x_train0=x[ind_train0,]
y_train0=y[ind_train0]

x_val0 = x[-ind_train0,]
y_val0 = y[-ind_train0]
t_val = t_all[-ind_train]
y_train0=y[ind_train0]


base_wf = workflow() %>% add_formula(y~.)
train_data_shift = as.tibble(x=x_train0)
train_data_shift[,'y'] = y_train0
base_shift_fit = base_wf %>%
  add_model(linear_reg() %>% set_engine("lm")) %>% fit(train_data_shift)  
cal_data_shift = as.tibble(x=x_val0)
cal_data_shift[,'y'] = y_val0
test_data_shift = as.tibble(x00)
test_data_shift[,'y'] = y00

quant_shift_base <-
  int_conformal_quantile(
    base_shift_fit, 
    train_data = train_data_shift,
    cal_data = cal_data_shift, 
    level = 0.90,
    ntree = 2000)

test_quant_shift_base <- 
  predict(quant_shift_base, test_data_shift) %>% 
  bind_cols(test_data_shift)
test_quant_shift_base
print(paste("lm: coverage", coverage(test_quant_shift_base), "width", round(width(test_quant_shift_base), digits=2)))
