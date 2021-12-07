# list of well-defined protocols for Data simulation

F.unif.int <- function(alpha){
  return(as.integer(runif(1,1,alpha)))
}

G.Bernoulli <- function(beta){
  return(rbinom(1,1,beta))
}

create.beta.matrix.bern <- function(K,p,a=1.0,b=1.0){
  if(length(a) == 1){
    a = rep(a,p)
  }
  else if(length(a) != p){
    stop("length of a is not either 1 or p")
  }
  if(length(b) == 1){
    b = rep(b,p)
  }
  else if(length(b) != p){
    stop("length of b is not either 1 or p")
  }
  return(sapply(1:p, function(i){rbeta(K,shape1=a[i],shape2=b[i])}))
}

#' DataSim
#'
#' A function that produces a matrix of fake single cell data following the
#' generative process outlined in the following figure.
#'
#' \if{html}{\figure{DataSim.png}{options: width=100 alt="Data sim fig"}}
#' \if{latex}{\figure{DataSim.png}{options: width=100 alt="Data Sim Figure"}}
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
#' TODO: finish the rest of this later
#'
#' @param x A number.
#' @param y A number.
#' @return
#' @examples
#'
#'
DataSim <- function(N,K,p,F,G,alpha,beta,seed=19890418,...){
  set.seed(seed)
  grp_cnts = rep(0,K)
  # if alpha is scalar with dim 1x1
  for(i in 1:(K-1)){
    grp_cnts[i] = min(F(alpha),N-sum(grp_cnts)-(K-i))
  }
  grp_cnts[K] = N-sum(grp_cnts)
  # if beta is kxp matrix
  data = matrix(0,nrow=N,ncol=p)
  grp_inclusion = rep(0,N)
  for(i in 1:K){
    for(j in 1:grp_cnts[i]){
      row_i = sum(grp_cnts[0:(i-1)])+j
      data[row_i,] = sapply(beta[i,], G)
      grp_inclusion[row_i] = i
    }
  }

  # store & return results
  result <- list("grp_cnts"=grp_cnts,"data"=data,"grp_inclusion"=grp_inclusion)
  return(result)
}


#' DataSim.unif.bern
#'
#' TODO: finish later
#'
#' @name DataSim.unif.bern
#'
#' @param N
#' The number of cells/rows to cluster
#' @return
#' @usage
#'
#'
DataSim.unif.bern <- function(N=100,K=10,p=2000,seed=19890418,beta.params=c(1.0,5.0),...){
  set.seed(seed)
  alpha <- 2*N/K
  beta <- create.beta.matrix.bern(K,p,beta.params[1],beta.params[2])
  return(DataSim(N,K,p,F=F.unif.int,G=G.Bernoulli,alpha=alpha,beta=beta,seed=seed))
}

