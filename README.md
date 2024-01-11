# Doubly-robust-prediction-under-covariate-shift

This repo contains the code and data for the simulation and real data application results in [this paper](https://arxiv.org/abs/2203.01761v3).  

Individual files do not have other dependencies with each folder containing a set of results, including the code to generate the results (mostly summary statistics such as coverage and size of prediction intervals), the .rda data results and the ggplot code to yield the plots.

For example, in ```/real_data_fig1```, ```real-500-absolute_0.9_0.1.R``` is used to produce the data ```real-500-absolute_0.9_0.1.rda```, which is then used by ```ggplot_real_all_0.9_0.9.R``` to produce the final figure ```ggplot_real_all_0.9_0.1.pdf```.