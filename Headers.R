## load R libraries we'll need for analysis 
library(reticulate)
library(tensorflow)
library(tfdatasets) ## format training dataset in tensorflow format
library(keras)      ## Deep learning packages, Keras API and tensorflow
library(tidyverse)  ## rbernoulli
library(data.table) ## 'shift' function
library(lubridate)
library(cvTools)
library(caret)
library(randomForestSRC) ## compute IPCW weights and IPCW pseudo RMST using random survival forest
library(pseudo)    ## for pseudosurv function in DNNSurv package
library(survival)

## plotting library
library(ggplot2)
library(cowplot)
library(gridExtra)
library(tidybayes)
