

# list of well-defined protocols for Data simulation

#' DataSim
#'
#' A function that produces a matrix of fake single cell data following the
#' generative process outlined in the following figure.
#'
#' \figure{DataSim.png}{options: width=100 alt="R logo"}
#' The user must specify values for N, p, and K, and either explicitly defines,
#' G()’s,F(),β’s, and ɑ or produces a generative process for them. First the
#' number of cells per group is decided by the function F() using parameters ɑ,
#' the last group has a fixed size to ensure the total number of cells is N.
#' Then for each cell the corresponding G(β) is used to produce the values of
#' the covariate assuming each cluster and gene has a different data generating
#' function G with different parameters β (Note: the user may choose to use the
#' same function; however, the model assumes they are different). The return
#' output would be a dataframe of N\*p size with a (k-1)\*1 vector detailing the
#' borders of the cell clusters by row.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#'
#' @export
DataSim <- function(N,K,p,F,G,alpha,beta,...){
  result <- list()
  return(result)
}
