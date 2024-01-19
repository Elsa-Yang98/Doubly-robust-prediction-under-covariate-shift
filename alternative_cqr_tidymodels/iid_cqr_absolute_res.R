# source('wcp_helper.R')
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
N=20000
p=4
beta = c(27.4,13.7,13.7,13.7)
dat.x=matrix(rnorm(N*p),N,p)
dat.y = 210+dat.x%*%beta[1:p]+rnorm(N)

train_x = dat.x[1:N/2,]; train_y = dat.y[1:N/2]
cal_x = dat.x[N/2:N*5/6, ]; cal_y = dat.y[N/2:N*5/6]
test_x = dat.x[N*5/6:N, ]; test_y = dat.y[N*5/6:N]

base_wf = workflow() %>% add_formula(y~.)
train_data = as.tibble(x=train_x)
train_data[,'y'] = train_y
base_fit = base_wf %>%
  add_model(linear_reg() %>% set_engine("lm")) %>% fit(train_data)  
cal_data = as.tibble(x=cal_x)
cal_data[,'y'] = cal_y
test_data = as.tibble(test_x)
test_data[,'y'] = test_y

quant_base <-
  int_conformal_quantile(
    base_fit, 
    train_data = train_data,
    cal_data = cal_data, 
    level = 0.90,
    ntree = 2000)

test_quant_base <- 
  predict(quant_base, test_data) %>% 
  bind_cols(test_data)
# test_quant_base
print(paste("cqr: coverage", coverage(test_quant_base), "width", round(width(test_quant_base), digits=2)))


split_int <- int_conformal_split(base_fit, cal_data)
# split_int
test_split_res <- 
  predict(split_int, test_data, level = 0.90) %>% 
  bind_cols(test_data)
print(paste("absolute residual: coverage", coverage(test_split_res), "width", round(width(test_split_res), digits=2)))

