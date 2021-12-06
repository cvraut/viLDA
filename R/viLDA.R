# viLDA.R
# this script wraps the stochastic variational inference c++ code

#' import stats
#' import Rcpp
#' import RcppArmadillo

#' viLDA.stoch
#'
#' TODO: fill out later
#'
#' @name viLDA.stoch
#' @param dat
#' The data in long format.
#' This is a data.frame with 3 columns:
#' $doc = doc ids
#' $words = word ids
#' $N = counts
#' @return
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
    v = length(unique(dat$words))
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
