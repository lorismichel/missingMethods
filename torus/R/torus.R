#' Creating an extra torus
#' @param X a matrix n \times p with possibly missing values encoded as NA
#' @param nb.nodes number of nodes in the torus.
#' TODO: add the possibility to choose a decision with NA or not
extraTorus <- function(X,
                       nb.nodes = 7) {

  if (nrow(na.omit(X))==0) {
    stop("cannot generate a torus.")
  }

  # creating the torus graph
  g <- igraph::graph.lattice(c(nb.nodes, nb.nodes),
                             directed = TRUE)

  # pasting the ends
  v1 <- seq(from = nb.nodes^2 - nb.nodes+1, to = nb.nodes^2, by = 1)
  v2 <- seq(from=1,to=nb.nodes,by=1)
  v3 <- seq(from = nb.nodes, to = nb.nodes^2, by = nb.nodes)
  v4 <- seq(from = 1, to = nb.nodes^2 - 1, by = nb.nodes)

  a  <- cbind(rbind(v1,v2), rbind(v3,v4))
  a2 <- matrix(a, nrow=length(a), ncol = 1)
  g  <- igraph::add.edges(g, a2)

  # get the adjacency matrix
  adj.mat <- as.matrix(igraph::as_adjacency_matrix(g))


  # randomize variables
  variable.mat <- adj.mat
  for (i in 1:nrow(adj.mat)) {
    variable.mat[i, adj.mat[i,]==1] <- sample(1:ncol(X),
                                              1,
                                              replace=TRUE)
  }

  # randomize splits
  split.mat <- adj.mat
  for (i in 1:nrow(adj.mat)) {
    v <- variable.mat[i,adj.mat[i,]==1][1]
    split.mat[i, adj.mat[i,]==1] <- sample(quantile(na.omit(X[,v]), probs = seq(0.001,0.999, length.out = 500)),1)#  sample(na.omit(X[,v]),1)
  }
 # split.mat[adj.mat==1] <- sapply(variable.mat[adj.mat==1], function(v) sample(na.omit(X[,v]),1))

  return(list(adj.mat = adj.mat, variable.mat = variable.mat, split.mat = split.mat))
}



#' @param X a matrix n \times p with possibly missing values encoded as NA
#' @param nb.nodes number of nodes in the torus.
#' TODO: add the possibility to choose a decision with NA or not
closedTree <- function(X,
                       depth = 2) {

  if (nrow(na.omit(X))==0) {
    stop("cannot generate a torus.")
  }

  # creating the tree graph
  g <- igraph::make_tree(n = 2^{depth+1}-1, mode = "out")

  # adj mat
  adj.mat <- as.matrix(igraph::as_adjacency_matrix(g))
  leaves <- which(apply(adj.mat,1,sum)==0)
  a <- rbind(leaves,1)
  a2 <- matrix(a, nrow=length(a), ncol = 1)
  g  <- igraph::add.edges(g, a2)

  # get the adjacency matrix after pasting the leaves
  adj.mat <- as.matrix(igraph::as_adjacency_matrix(g))


  # which are not leaves
  nb.no.leaves <- length(which(apply(adj.mat,1,sum)!=1))

  # randomize variables
  variable.mat <- adj.mat
  for (i in 1:nb.no.leaves) {
    variable.mat[i, adj.mat[i,]==1] <- sample(1:ncol(X),
                                              1,
                                              replace=TRUE)
  }

  # randomize splits
  split.mat <- adj.mat
  for (i in 1:nb.no.leaves) {
    v <- variable.mat[i,adj.mat[i,]==1][1]

    split.mat[i, adj.mat[i,]==1] <- sample(quantile(na.omit(X[,v]), probs = seq(0.01,0.99, length.out = 1000)),1)

  }
  #split.mat[adj.mat==1] <- sapply(variable.mat[adj.mat==1], function(v) sample(na.omit(X[,v]),1))

  return(list(adj.mat = adj.mat, variable.mat = variable.mat, split.mat = split.mat))
}

#' Get the stationary distribution over a torus
#' @param X a matrix of observatiosn, n \times p (with NA possibly)
#' @param method method to compute stationary distribution, either "power" or "eigen"
#' @param prob a 3-dim numeric vector with (p_s, p_r, 1-p_s-p_r) where p_s is probability of staying, p_r choosing random neighbor and the rest for data-driven decision
#' @param power the power used in method "power"
stationaryDistr <- function(object,
                            X,
                            method = "power",
                            prob = c(0.05, 0.05, 1-0.05-0.05),
                            power = 100,
                            subset = NULL) {

  # if we have missing values, we use NA
  X[!is.finite(X)] <- NA

  pi.mat <- t(apply(X, 1, function(x) {

    tr.mat <- getTransitionMatrix(object,
                                  x,
                                  prob = prob,
                                  power = power)

    # getting the stationary distribution
    if (method == "eigen") {
        e <- eigen(t(tr.mat))
        id <- which.min(abs(e$values-1))
        v <- e$vectors[,id]
        pi.stat <- abs(v) / sum(abs(v))
      } else {
        tr.mat.power <- powerMat(tr.mat, power = power)

        if (any(apply(tr.mat.power, 2, sd)>=0.001)) {
          warning("use method eigen.")
        }
        pi.stat <- apply(tr.mat.power, 2, mean)
      }

      return(pi.stat)
  }))

  if (!is.null(subset)) {
    pi.mat <- pi.mat[,subset]
    pi.mat <- pi.mat / rowSums(pi.mat)
  }

  return(pi.mat)
}


getTransitionMatrix <- function(object,
                                x,
                                prob = c(0.05, 0.05, 1-0.05-0.05),
                                power = 100) {

  # if we have missing values, we use NA
  x[!is.finite(x)] <- NA

  tr.mat <- object$adj.mat

  # nb no leaves
  nb.no.leaves <- length(which(apply(tr.mat, 1, sum)!=1))
    # data-driven transition + random allocation
  for (i in 1:nb.no.leaves) {

      # get the x value at the variable of the node
      x.val <- x[object$variable.mat[i, object$adj.mat[i,]==1][1]]
      # get the decision direction
      cond <- x.val <= object$split.mat[i,object$adj.mat[i,]==1][1]

      tr.mat[i, object$adj.mat[i,]==1][1] <- prob[3]*ifelse(is.na(x.val), 1/2, as.numeric(cond)) + prob[2] * 1/2
      tr.mat[i, object$adj.mat[i,]==1][2] <- prob[3]*ifelse(is.na(x.val), 1/2, 1-as.numeric(cond)) + prob[2] * 1/2
    }

    if (nb.no.leaves!=nrow(tr.mat)) {
      for (i in (nb.no.leaves+1):nrow(tr.mat)) {
        tr.mat[i,] <- tr.mat[i,]*(prob[2]+prob[3])
      }
    }

    # adding the "staying" transition
    tr.mat <- tr.mat + prob[1]*diag(1,nrow(tr.mat),ncol(tr.mat))

    return(tr.mat)
}


getSampleWeightsEnsemble <- function(object,
                             X,
                             d="inner_product",
                             subset = NULL,
                             method = "power",
                             prob = c(0.05, 0.05, 1-0.05-0.05),
                             power = 100) {

  l <- lapply(object, function(o) {
    # getting the stationary distr.
    pi.mat <- stationaryDistr(object = o, X = X,method=method,prob=prob, power=power,subset = subset)

    # computing norm
    #w <- apply(pi.mat, 1, function(x1) apply(pi.mat, 1, function(x2) sqrt(sum((x1-x2)^2)) ))

    # computing the kernel
    #w <- apply(pi.mat, 1, function(x1) apply(pi.mat, 1, function(x2) sum(x1*x2) )) #
    # normalizing the kernel to have empirical distributions
    w <- distance(pi.mat , method=d)

    return(w)
  })

  w <- Reduce(l, f = function(x,y) x+y)
  w <- w / length(l)
  w <- w / rowSums(w)
}


getSampleWeights <- function(object,
                             X,
                              d="inner_product" ,
                             subset = NULL,
                             ...) {
    # getting the stationary distr.
    pi.mat <- stationaryDistr(object = object, X = X, subset = subset,...)


    #w <- apply(pi.mat, 1, function(x1) apply(pi.mat, 1, function(x2) sum(x1*x2)))
    w<-distance(pi.mat , method=d)


    # normalizing the kernel to have empirical distributions

    #if (d=="inner_product"){
    w <- w / rowSums(w)
    #}

    return(w)

}


combine <- function(...) {
  return(list(...))
}


