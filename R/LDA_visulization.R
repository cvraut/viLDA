#'LDA_visualization

#'@param X: the matrix include doc word and N
#'@param Cluster_Label: Predicted Label
#'@param Method: normalized method
#'@param dims: dimension to visualized
#'@param Three_dimension_Plot: if dims = 3, whether generate a 3D plot or scatterplot for each pair of dimensions.
#'@param seed:
#'@return visualization plot
#'
#'@import ggpubr,ggplot2, Rtsne, rgl
#'@export



LDA_visulization <- function(X,
                             Cluster_Label,
                             Method = c("proportion","mean","median"),
                             dims = c(2,3),
                             Three_dimension_Plot = FALSE,
                             seed = 12345){
  # ensure the dims are 2 or 3
  stopifnot("Dims should be either 2 or 3" = dims %in% c(2,3))

  # transform to the count matirx
  Count <- X %>%
    pivot_wider(names_from = word,
                values_from = N) %>%
    select(-doc)
  Count <- as.matrix(Count)
  Count <- ifelse(is.na(Count),0,Count)
  Count <- t(Count)


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
    res <- .vi_lda(Count,Method,dims)
  }else{
    set.seed(seed)
    res <- .vi_lda(Count,Method,dims)
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
    if(Three_dimension_Plot){
      plot3d(
        x=final$V1, y=final$V2, z=final$V3,
        col = final$Cluster,
        type = 's',
        radius = 1,
        xlab="Dimension 1", ylab="Dimension 2", zlab="Dimension 3")
    }else{
      p1 = ggplot(final,aes(V1,V2,col = Cluster)) + geom_point() +
        labs(x = "Dimension 1",
             y = "Dimension 2")
      p2 = ggplot(final,aes(V1,V3,col = Cluster)) + geom_point() +
        labs(x = "Dimension 1",
             y = "Dimension 3")
      p3 = ggplot(final,aes(V2,V3,col = Cluster)) + geom_point() +
        labs(x = "Dimension 2",
             y = "Dimension 3")
      list(p1,p2,p3)
    }
  }
}
