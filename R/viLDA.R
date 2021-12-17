# viLDA.R
# this script wraps the stochastic variational inference c++ code

#' viLDA.stoch
#'
#' @param dat
#' The data in long format.
#' This is a data.frame with 3 columns:
#' $doc = doc ids
#' $word = word ids
#' $N = counts
#' An example of this format can be generated through viLDA::data_gen
#'
#' @param K
#' Integer, the number of clusters
#'
#' @param d
#' Intger, the number of documents (default = NULL)
#' if NULL then this is imputed from the dat object
#'
#' @param v
#' Integer, the size of the vocabulary (default = NULL)
#' if NULL then this is also imputed from the dat object
#'
#' @param maxIterVal
#' Integer, the maximum number of outer iterations for the global optima search (default = NULL)
#' if NULL then this imputed to be log(v)
#'
#' @param maxVBIterVal
#' Integer, the maximum number of inner iterations for the local optima search (default = NULL)
#' if NULL then this is imputed to be ceiling(sqrt(v))
#'
#' @param alphaWords
#' Float, the prior for the word-cluster distributions (default = 0.2)
#'
#' @param alphaTopics
#' Float, the prior for the document-cluster distributions (default = 0.2)
#'
#' @param rho
#' Float, the learning rate for the parameter updates (default = 0.001)
#'
#' @param tol
#' Float, the tolerance for assessing convergence (default = 0.1)
#'
#' @param seed
#' Float, the seed run at the begining (default = 19890418)
#'
#' @param ...
#' extra arguments can be passed, but they don't do anything yet 🛠
#'
#' @return list([[1]],[[2]],[[3]])
#'  - [[1]]: k x v matrix of the expected values for the word-cluster probabilities
#'  - [[2]]: I x k matrix of the complete stored values of the parameters through all iterations
#'  - [[3]]: n length vector of the true value
#'
#' @export
viLDA.stoch <- function(dat,
         K,
         d=NULL,
         v=NULL,
         maxIterVal=NULL,
         maxVBiterVal=NULL,
         alphaWords = 0.2,
         alphaTopics=0.2,
         rho = 0.001,
         tol=0.1,
         seed=19890419,...){
  set.seed(seed)
  if(is.null(d)){
    d = length(unique(dat$doc))
  }
  if(is.null(v)){
    v = length(unique(dat$word))
  }
  if(is.null(maxIterVal)){
    maxIterVal <- log(v)
  }
  if(is.null(maxVBiterVal)){
    maxVBiterVal <- ceiling(sqrt(log(v)))
  }

  resList <- svi(data=as.matrix(dat),
                 numDistinctWordVec=as.vector(table(dat$doc)),
                 topics = 0:(K-1),
                 lengthVocab = v,
                 numDocuments = d,
                 maxIterConst = maxIterVal,
                 maxVBiterConst = maxVBiterVal,
                 alphaWords = alphaWords,
                 alphaTopics = alphaTopics,
                 rho = rho,
                 tol = tol)
  return(resList)
}
