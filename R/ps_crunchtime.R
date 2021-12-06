source("R/data_gen.R")
source("R/assessments.R")

# libs
library(stats)
library(Rcpp)
library(RcppArmadillo)
library(tictoc)
library(topicmodels)
library(tidyverse)
library(tidytext)

set.seed(56)
d = 1024
v = 128
k = 2
data_obj <- data_gen(d,v,k, doc_length_scale = 32, voc_p_scale = 8, spike_overlap = 0.05)
# current best
dat_matrix = data.frame(data_obj$dat) %>% cast_dtm(doc,word,N)
tic()
mod = topicmodels::LDA(dat_matrix,k = 2, method = "VEM")
tm_res = posterior(mod,dat_matrix)
toc()
true = data_obj$gen_topics+1
pred = get_plurarity_topics(tm_res$topics)
get_best_mapping(pred,true)$prop_correct

# ours
sourceCpp("viAlgorithm.cpp")
tic()

### VERY SMALL ITERATION VALUES BUT IT STILL WORKS
maxIterVal <- log(v)
maxVBiterVal <- ceiling(sqrt(log(v)))

### USE SMALL LEARNING RATE (rho)
### MAKE ITERATION VALUES VERY SMALL (still manages to converge fast)
resList <- svi(data=as.matrix(data_obj$dat),
               numDistinctWordVec=as.vector(table(data_obj$dat$doc)),
               topics = 0:(k-1),
               lengthVocab = v,
               numDocuments = d,
               maxIterConst = maxIterVal,
               maxVBiterConst = maxVBiterVal,
               alphaWords = 0.2,
               alphaTopics = 0.2,
               rho = 0.001,
               tol = 0.1)
toc()
true = data_obj$gen_topics+1
pred = get_plurarity_topics(resList[[3]])
print(get_best_mapping(pred,true)$prop_correct)


