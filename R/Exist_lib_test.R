library(topicmodels)
library(tidyverse)
library(tidytext)
library(tictoc)
library(doParallel)
library(gtools)


location <- getwd()
d <- c(500,1000,2000)
v <- c(100,500,1000)
t <- c(2,5,10)
para_dataset <- expand.grid(d,v,t)
myCluster <- makeCluster(6,
                         type = "PSOCK")

registerDoParallel(myCluster)

topic_model_time <- foreach(i = 1:3^3, .combine = rbind) %dopar% {
  library(topicmodels)
  library(tidyverse)
  library(tidytext)
  library(tictoc)
  temp <- load(paste(location,
             "/test_data/data_file_d",
             para_dataset[i,1],
             "_v",
             para_dataset[i,2],
             "_t",
             para_dataset[i,3],
             ".Rda",sep = ""))

  tic()
  lda <- topicmodels::LDA(data_obj[[1]] %>%
                            cast_dtm(doc, word,N),
                          control = list(alpha = 0.2, iter = 2000),
                          k = para_dataset[i,3],
                          method= "Gibbs")
  t1 <- toc()
  print(c(para_dataset[i,],as.numeric(t1$toc - t1$tic)))
}

topic_model_time <- data.frame(topic_model_time)
colnames(topic_model_time) <- c("d","v","t","time")
# save(topic_model_time,file = "topic_model_time.rdata")
load("topic_model_time.rdata")



LDA_model_time <- foreach(i = 1:3^3, .combine = rbind) %dopar% {
  library(lda)
  library(tidyverse)
  library(tictoc)
  load(paste(location,
              "/test_data/data_file_d",
               para_dataset[i,1],
              "_v",
              para_dataset[i,2],
              "_t",
              para_dataset[i,3],
              ".Rda",sep = ""))
  tic()
  x <- sapply(0:(para_dataset[i,1]-1),FUN = function(i){
    temp = data_obj[[1]]
    temp = temp[temp[,1] == i,2:3] %>%
      arrange(word)
    rownames(temp) = temp[,1] - 1
    colnames(temp) = NULL
    t(temp)
  },
  simplify = FALSE)
  vocab_vec <- 0:(para_dataset[i,2]-1)
  t1 <- toc()


  tic()
  lda_lda <- lda::lda.collapsed.gibbs.sampler(x,
                                        vocab_vec,
                                        num.iterations = 2000,
                                        K = para_dataset[i,3],
                                        alpha = 0.2,
                                        eta = 0.2)
  t2 <- toc()
  print(c(para_dataset[i,],
            as.numeric(t1$toc - t1$tic),
            as.numeric(t2$toc - t2$tic)))

}

LDA_model_time <- data.frame(LDA_model_time)
colnames(LDA_model_time) <- c("d","v","t","time")
# save(LDA_model_time,file = "LDA_model_time.rdata")
load("topic_model_time.rdata")


#-----------------------------------#
# test for svi#

library(stats)
library(Rcpp)
library(RcppArmadillo)
library(tictoc)
library(topicmodels)
library(tidyverse)
library(tidytext)


i = 27 # for 2000-1000-10
load(paste(location,
           "/test_data/data_file_d",
           para_dataset[i,1],
           "_v",
           para_dataset[i,2],
           "_t",
           para_dataset[i,3],
           ".Rda",sep = ""))

sourceCpp("viAlgorithm.cpp")


resList <- svi(data = as.matrix(data_obj[[1]]),
               numDistinctWordVec= as.vector(table(data_obj[[1]][,1])),
               topics = 0:9,
               lengthVocab = 1000,
               numDocuments = 2000,
               maxIterConst = ceiling(log(1000)),
               maxVBiterConst = ceiling(log(1000)),
               alphaWords = 0.2,
               alphaTopics = 0.2,
               rho = 0.1,
               tol = 0.0001)
pred = get_plurarity_topics(resList[[3]])





tic()
result_map <- get_best_mapping(pred,data_obj$gen_topics + 1)
toc()


result_map$prop_correct



