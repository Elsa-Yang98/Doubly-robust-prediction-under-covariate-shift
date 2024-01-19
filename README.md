# Doubly-robust-prediction-under-covariate-shift

This repo contains the code and data for the simulation and real data application results in [this paper](https://arxiv.org/abs/2203.01761v3).  

- `conformal.pred.sl.fast()` is the function to implement Weighted Conformal Prediction (WCP), which is modified upon the code available at [http://www.github.com/ryantibs/conformal/](http://www.github.com/ryantibs/conformal/) for the original [WCP paper](https://arxiv.org/abs/1904.06019)


Individual files do not have other dependencies with each folder containing a set of results, including the code to generate the results (mostly summary statistics such as coverage and size of prediction intervals contained in ```.rda``` files). The files starting with ```ggplot_``` are used to yield the plots illustrated in the paper. 

##
For example, in ```/real_data_fig1```, ```real-500-absolute_0.9_0.1.R``` is used to produce the data ```real-500-absolute_0.9_0.1.rda```, which is then used by ```ggplot_real_all_0.9_0.9.R``` to produce the final figure ```ggplot_real_all_0.9_0.1.pdf```.