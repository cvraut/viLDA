# assessments.R
# Code to assess performance and evaluate results of viLDA
#' import gtools

# hidden:
# maps predicted values according to the specified mapping
map_pred <- function(mapping,pred){
  mapped_predict = rep(NA,length(predict))
  for(i in 1:length(mapping)){
    j = mapping[i]
    mapped_predict[pred == i] = j
  }
  return(mapped_predict)
}
# returns the proportion of correctly classified values according to
# predictions and true mappings
prop_correctly_classified <- function(mapping,predict,true_val){
  mapped_predict = map_pred(mapping,predict)
  correct = 0
  for(i in 1:length(true_val)){
    if(true_val[i]==mapped_predict[i]){
      correct = correct+1
    }
  }
  return(correct/length(true_val))
}

# unhidden:
#' get_plurarity_topics
#'
#' get the 1-indexed topics by plurarity of the document to topic matrix
#' @name get_plurarity_topics
#'
#' @param doc_2_top_mat
#' n by k matrix of floating point numbers
#' The rows should sum to 1
#'
#' @return n-length vector of plurarity mappings
#' The values will be integers ranging from [1-K]
#'
#' @usage get_plurarity_topics(doc_2_top_mat)
#'
#' @export
get_plurarity_topics <- function(doc_2_top_mat){
  return(apply(doc_2_top_mat,1,FUN=which.max))
}

#
# TODO: finish documenting later :/
#' get_best_mapping
#'
#' returns the best mapping through complete search
#' @name get_best_mapping
#'
#' @param predicted
#' n length vector of predicted groups
#'
#' @param true
#' n length vector of true groups
#'
#' @return list("prop_correct","mapped_pred","true_val","mapping")
#' $prop_correct: float [0-1] proportion of the true value, that best prediction map correctly gets
#' $mapped_pred: n length vector, mapped predictions of best mapping
#' $true_val: n length vector of the true value
#' $mapping: k length vector (k is the number of groups) of the best group mappings
#'
#' @usage get_best_mapping(predicted,true)
#'
#' @export
get_best_mapping <- function(predicted,true){
  K = max(true)
  all_mappings = permutations(n = K, r = K, v = 1:K)
  results = apply(all_mappings,1,FUN = function(perm){prop_correctly_classified(perm,predicted,true)})
  max_i = which.max(results)
  return(list("prop_correct"=results[max_i],"mapped_pred"=map_pred(all_mappings[max_i,],predicted),"true_val"=true,"mapping"=all_mappings[max_i,]))
}
