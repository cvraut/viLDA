#' import stats
#' import Rcpp
#' import RcppArmadillo

#' viLDA.gibbs
#'
#' TODO: fill out later
#'
#' @name viLDA.gibbs
#' @param N
#' The number of documents
#'
#' @return
#'
#' @export
viLDA.gibbs <- function(N,K,p,data,alpha,beta,seed=19890419,...){
  set.seed(seed)
  data.long = c(t(data))
  data.id = c(sapply(1:N,FUN=function(i){rep(i - 1, p)}))
  posteriorSamples <- gibbsSampler(words = data.long,
                                   docIDs = docIDs,
                                   topics = 0:(K-1),
                                   lengthVocab = p,
                                   numDocuments = N,
                                   alphaWords = beta,
                                   alphaTopics = alpha,
                                   numEpochs = 50000,
                                   warmUp = 5000,
                                   lag = 50)

  # Obtain posterior means
  posteriorAlphaWords <- Reduce("+", posteriorSamples[[1]])/length(posteriorSamples[[1]])
  posteriorAlphaTopics <- Reduce("+", posteriorSamples[[2]])/length(posteriorSamples[[2]])
  result <- list("posterior_gene_to_cluster"=posteriorAlphaWords,
                 "posterior_cell_to_cluster"=posteriorAlphaTopics,
                 "posterior_samples"=posteriorSamples)
  return(result)
}
