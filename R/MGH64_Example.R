library(topicmodels)
library(tidyverse)
library(tidytext)
library(tictoc)
library(Rcpp)
library(RcppArmadillo)
source("R/LDA_visulization.R")
sourceCpp("viAlgorithm.cpp")
source("R/assessments.R")
load("data/MGH64.RData")

#-------------------------# remove 0
MGH64_update = MGH64[rowSums(MGH64) != 0,]
MGH64_update = MGH64_update[,colSums(MGH64_update) != 0]
mat_train = as.matrix(t(MGH64_update))
sum(rowSums(mat_train) == 0)
sum(colSums(mat_train) == 0)

document_matrix <- data.frame(t(mat_train))
document_matrix$document <- rownames(document_matrix)
document_matrix <- document_matrix %>%
  pivot_longer(cols = !document,
               names_to = "term",
               values_to = "count") %>%
  filter(count > 0)
colnames(document_matrix) = c("doc","word","N")

topic_model_matrix <- document_matrix  %>%
  cast_dtm(doc, word, N)





tic()
lda_topic <- topicmodels::LDA(topic_model_matrix, control = list(alpha = 0.2, iter = 2000), k = 2, method= "Gibbs")
toc()

# save(lda_topic,file = "topic_lda_MGH64.rda")
pred_lda <- posterior(lda_topic,topic_model_matrix)
pred_lda_label <- apply(pred_lda$topics,1,FUN = function(i) which.max(i))


LDA_visulization(document_matrix,pred_lda_label,Method = c("proportion"),seed = 12345)



training_matrix_generator <- function(X = NULL, Words = NULL, Documents = NULL, count_matrix = TRUE){
  Words = Words
  Documents = Documents
  if(count_matrix){
    # each row of count_matrix represent the number of documents
    # each col of count_matrix represent the number of words
    X = data.frame(X)
    Unique_docu_name = rownames(X)
    Unique_word_name = colnames(X)
    if(is.null(Unique_docu_name)){
      Unique_docu_name = 1:nrow(X)
    }
    if(is.null(Unique_word_name)){
      Unique_word_name = 1:ncol(X)
    }
    X$Documents = Unique_docu_name
    X <- X %>%
      pivot_longer(cols = !Documents,
                   names_to = "Words",
                   values_to = "n") %>%
      filter(n!=0)
  }else{
    X = data.frame(Documents = Documents, Words = Words)
    X = X %>%
      group_by(Documents) %>%
      count(Words) %>%
      ungroup() %>%
      filter(n != 0)
    Unique_word_name = sort(unique(Words))
    Unique_docu_name = sort(unique(Documents))
  }
  word_index = 1:length(Unique_word_name) - 1
  Docu_index = 1:length(Unique_docu_name) - 1
  X$Documents = match(X$Documents, Unique_docu_name) - 1
  X$Words = match(X$Words,Unique_word_name) - 1
  return(list(mapping_docu = data.frame(documents = Unique_docu_name,
                                        index = Docu_index),
              mapping_word = data.frame(words = Unique_word_name,
                                        index = word_index),
              X = X))

}

svi_matrix <- training_matrix_generator(X = t(mat_train))

tic()
resList <- svi(data = as.matrix(svi_matrix$X),
               numDistinctWordVec= as.vector(table(svi_matrix$X[,1])),
               topics = 0:1,
               lengthVocab = nrow(unique(svi_matrix$X[,2])),
               numDocuments = nrow(unique(svi_matrix$X[,1])),
               maxIterConst =  ceiling(log(nrow(unique(svi_matrix$X[,2])))),
               maxVBiterConst = ceiling(log(nrow(unique(svi_matrix$X[,2])))),
               alphaWords = 0.2,
               alphaTopics = 0.2,
               rho = 1/100000,
               tol = 0.000000001)
toc()

pred_svi = get_plurarity_topics(resList[[3]])
colnames(svi_matrix$X) = c("doc","word","N")

LDA_visulization(document_matrix,pred_svi,dims = 2,Method = c("proportion"),seed = 12345)
