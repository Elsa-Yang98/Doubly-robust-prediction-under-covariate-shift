##########################code from site##########################
library(tidymodels)
tidymodels_prefer()

set.seed(1)
make_data <- function(n, std_dev = 1 / 5) {
  tibble(x = runif(n, min = -1)) %>%
    mutate(
      y = (x^3) + 2 * exp(-6 * (x - 0.3)^2),
      y = y + rnorm(n, sd = std_dev)
    )
}

n <- 1000
set.seed(8383)
train_data <- make_data(n)
cal_data  <- make_data(250)
test_data <- make_data(10000)            


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



set.seed(1)
base_wf = workflow() %>%
  add_formula(y~x)
base_fit = base_wf %>%
  add_model(linear_reg() %>%
              set_engine("lm")) %>%
  fit(train_data)

quant_int_base <-
  int_conformal_quantile(
    base_fit, 
    train_data = train_data,
    cal_data = cal_data, 
    level = 0.90,
    ntree = 2000)
test_quant_base <- 
  predict(quant_int_base, test_data) %>% 
  bind_cols(test_data)
test_quant_base
print(paste("lm: coverage", coverage(test_quant_base), "width", round(width(test_quant_base), digits=2)))

##########################my cqr##########################
set.seed(1)
quantiles = c(0.05, 0.95)
fit_r = grf::quantile_forest(as.matrix(train_data[,'x']),as.matrix(train_data[,'y']), quantiles = quantiles)
# quantile_train0 = predict(fit_r, train_data[,'x'], quantiles = quantiles)$predictions
# r_train0 = as.matrix(pmax( quantile_train0[,1]- train_data[,'y'],train_data[,'y']- quantile_train0[,2] ))
quantile_val = predict(fit_r, cal_data[,'x'], quantiles = quantiles)$predictions
r_val = as.matrix(pmax( quantile_val[,1]- cal_data[,'y'], cal_data[,'y']- quantile_val[,2] ))

wts = rep(1,nrow(train_data[,'x'])+nrow(cal_data[,'x'])+nrow(test_data[,'x']))
pred_r_new = predict(fit_r, data.frame(test_data[,'x']), quantiles = quantiles)$predictions
out = conformal.pred.sl.fast(rbind(as.matrix(train_data[,'x']),as.matrix(cal_data[,'x'])),as.matrix(test_data[,'x']), ind_train0=1:nrow(train_data[,'x']),m=1,pred=pred_r_new, res=r_val,w=wts, type ="CQR", quantiles = quantiles)
wcp_cov_lm = mean(out$lo <= test_data[,'y'] & test_data[, 'y']<= out$up)
wcp_len_lm = mean(out$up - out$lo)



