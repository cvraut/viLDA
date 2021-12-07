# data_gen.R
# scripts to generate test data for viLDA

#' data_gen
#'
#' Technically this returns the number of "code points", in a string. One
#' code point usually corresponds to one character, but not always. For example,
#' an u with a umlaut might be represented as a single character or as the
#' combination a u and an umlaut.
#'
#' @name data_gen
#' @param n_doc
#' The number of documents to generate
#' @return
#' @usage
#'
#' @importFrom data.table as.data.table
#' @export
data_gen <- function(n_doc,
                     n_vocab,
                     n_top,
                     doc_length_scale = 8,
                     doc_length_scale_var = 2,
                     voc_p_scale = 4,
                     spike_overlap = 0.05,
                     alphaWords = 0.2,
                     alphaTopics = 0.2,
                     seed=19890418,
                     topic_mix=FALSE,
                     DEBUG=FALSE){
  set.seed(seed)
  vocab = seq(0,n_vocab-1)
  lengthDocuments <- 1+rbinom(n_doc,size=doc_length_scale*n_vocab,prob=1/doc_length_scale_var)

  numWords <- sum(lengthDocuments)
  topics <- 0:(n_top-1)

  # initialize word distributions
  z = matrix(0,nrow=n_top,ncol=n_vocab)
  z0 <- rep(1/(voc_p_scale*n_vocab),n_vocab)
  centers <- rep(NA,n_top)
  centers[1] <- 3
  center_shift = -2*qnorm(spike_overlap/2)
  for(i in 2:n_top){
    centers[i]<- centers[i-1]+center_shift
  }
  norm_length <- centers[n_top]+3
  norm_scale <- n_vocab/norm_length

  for(k in 1:n_top){
    z0 <- rep(1/(voc_p_scale*n_vocab),n_vocab)
    z0_help <- rep(0,n_vocab)
    z0_help[1] <- pnorm(1/norm_scale-(3*k))
    for(rhs in 2:n_vocab){
      z0_help[rhs] <- (pnorm(rhs/norm_scale-centers[k])-pnorm((rhs-1)/norm_scale-centers[k]))
    }
    z0_help <- z0_help/(sum(z0_help))*(1-sum(z0))
    z0 = z0+z0_help
    z[k,] = z0
  }
  if(DEBUG){
    print(z)
  }
  wordDistributions <- z

  if(topic_mix){
    generatedTopics <- lapply(1:n_doc, function(x) x=c())
  } else {
    generatedTopics <- rep(NA, n_doc)
  }
  doc <- c()
  word <- c()
  count <- c()

  for (i in 1:n_doc) {
    doc <- c(doc, rep(i - 1, lengthDocuments[i]))
    if(topic_mix){
      z = sample(topics, lengthDocuments[i])
      sampledWords <- sapply(1:lengthDocuments[i],FUN=function(j){sample(vocab, 1, replace = TRUE, prob = wordDistributions[z[j]+1,])})
    } else {
      z = sample(topics, 1)
      sampledWords <- sample(vocab, lengthDocuments[i], replace = TRUE, prob = wordDistributions[z+1,])
    }
    generatedTopics[i] <- z
    word <- c(word, sampledWords)
  }
  dat = data.frame("doc"=doc,"word"=word)
  dat = dat[order(dat$doc, dat$word),]
  dat = data.frame(as.data.table(dat)[, .N, by = c('doc','word')])
  return(list("dat"=dat,"word_dist"=wordDistributions,"gen_topics"=generatedTopics,"doc_len"=lengthDocuments))
}

