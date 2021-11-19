library(celda)
library(SingleCellExperiment)


# Produce training data
# -----------------------------------------------

load(file="MGH64.RData")

#----------------------------------------------
# remove cell and feature which have all 0 -- row sum or col sum should be greater than 0
MGH64_update = MGH64[rowSums(MGH64) != 0,]
MGH64_update = MGH64_update[,colSums(MGH64_update) != 0]
mat_train = as.matrix(t(MGH64_update))
sum(rowSums(mat_train) == 0)
sum(colSums(mat_train) == 0)

#----------------------------#

# example
MGH_K_5 <- celda_C(x = mat_train,
                   K = 5,
                   nchains = 5,
                   algorithm = "EM",
                   zInitialize = "random",
                   splitOnIter = 1000)


altExp_featrue <- altExp(MGH_K_5,"featureSubset")
Cluster_Label <- altExp_featrue$celda_cell_cluster
#-----------------------------#



library(Rtsne)
library(ggplot2)
library(ggpubr)
# @count: count_matrix with dimension n * p
# @Cluster_Label: Label matrix with n * 1
# @method: normalization method
# @dims: dimension of vector
LDA_visulization <- function(count,
                             Cluster_Label,
                             method = c("proportion","mean","median"),
                             dims = c(2,3),
                             seed = 12345){
  .vi_lda <- function(count,
                      method,
                      dims){
    final <- matrix(NA, nrow = ncol(count), ncol = dims)
    cs <- colSums(count)
    norm <- switch(
      method,
      "proportion" = sweep(count, 2, cs, "/"),
      "median" = sweep(count, 2, cs / stats::median(cs), "/"),
      "mean" = sweep(count, 2, cs / mean(cs), "/")
    )
    res <- Rtsne::Rtsne(
      t(norm),
      dims = dims,
      pca = TRUE,
      max_iter = 2500,
      perplexity = 20,
      check_duplicates = FALSE,
      is_distance = FALSE,
      theta = 0,
      initial_dims = 20)
    return(res)
  }

  if(is.null(seed)){
      res <- .vi_lda(count,method,dims)
  }else{
    set.seed(seed)
    res <- .vi_lda(count,method,dims)
  }

  final <- res$Y
  if(dims == 2){
    Cluster_Label <- as.factor(as.character(Cluster_Label))
    final <- data.frame(V1 = final[,1],
                        V2 = final[,2],
                        Cluster = Cluster_Label)
    ggplot(final,aes(V1,V2,col = Cluster)) + geom_point() +
      labs(x = "Dimension 1",
           y = "Dimension 2")
  }else{
    Cluster_Label <- as.factor(as.character(Cluster_Label))
    final <- data.frame(V1 = final[,1],
                        V2 = final[,2],
                        V3 = final[,3],
                        Cluster = Cluster_Label)
    p1 = ggplot(final,aes(V1,V2,col = Cluster)) + geom_point() +
      labs(x = "Dimension 1",
           y = "Dimension 2")
    p2 = ggplot(final,aes(V1,V3,col = Cluster)) + geom_point() +
      labs(x = "Dimension 1",
           y = "Dimension 3")
    p3 = ggplot(final,aes(V2,V3,col = Cluster)) + geom_point() +
      labs(x = "Dimension 2",
           y = "Dimension 3")
    ggarrange(p1,p2,p3,nrow = 3)
  }
}

LDA_visulization(mat_train,
                 Cluster_Label,
                 method = "proportion",
                 dims = 3)

