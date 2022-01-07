#install.packages("devtools")
library(devtools)
library(roxygen2)

my.Rpackage <- as.package("../viLDA")
document(my.Rpackage)
load_all(my.Rpackage)
build()
install()

#devtools::install_github("cvraut/viLDA")

library(viLDA)
library(tictoc)

set.seed(56)
d = 1024
v = 128
k = 2
data_obj <- data_gen(d,v,k, doc_length_scale = 32, voc_p_scale = 8, spike_overlap = 0.05)

tic()
resList <- viLDA.stoch(dat=data_obj$dat,
                       K=2,
                       alphaWords = 0.2,
                       alphaTopics = 0.2,
                       rho = 0.001,
                       tol = 0.1)
true = data_obj$gen_topics+1
toc()
pred = get_plurarity_topics(resList[[3]])
print(get_best_mapping(pred,true)$prop_correct)
