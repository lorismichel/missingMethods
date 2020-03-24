# Random node graph 

# architecture of graph nodes
#     7    8    9
#     |    |    -
#     -    -    |
#3 -> 1 -> 2 -> 3 -> 1
#     |    |    -
#     -    -    |
#6 <- 4 -> 5 -> 6 <- 4
#     -    |    |
#     |    -    -
#9 -> 7 <- 8 <- 9 -> 7
#     |    |    -
#     -    -    |
#     1    2    3



############################ torus graph

# igraph version
require("igraph")
require("rgl")
n <- 7
g <- graph.lattice(c(n,n), directed = TRUE) # create a square lattice (nxn)

plot(g,vertex.size = 0.5,vertex.size = 4,vertex.label = 1:9,vertex.color = "red")
# want to connect up the corners (horribly done)
v1 <- seq(from = n^2 - n+1, to = n^2, by = 1)
v2 <- seq(from=1,to=n,by=1)
v3 <- seq(from = n, to = n^2, by = n)
v4 <- seq(from = 1, to = n^2 - 1, by = n)

a <- cbind(rbind(v1,v2), rbind(v3,v4))
a2 <- matrix(a,nrow=length(a),ncol = 1)
g <- add.edges(g,a2)
plot(g,vertex.size = 0.5,vertex.size = 4,vertex.color = "red",
     vertex.label = 1:9)

graph.mat <- as.matrix(as_adjacency_matrix(g))
rowSums(graph.mat)
colSums(graph.mat)


############################ data

# data
n <- 200 # number of observations
d <- 5 # number of dimensions

X <- MASS::mvrnorm(n = n, mu = rep(0, d), 
                   Sigma = diag(0.6,d,d)+matrix(0.4,d,d)) # naive matrix example
Xtest <- MASS::mvrnorm(n = n, mu = rep(0, d), 
                       Sigma = diag(0.6,d,d)+matrix(0.4,d,d)) # naive matrix example

# NA filter
na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(0.95, 0.05), replace = TRUE))
Xna <- X
Xna[na==1] <- NA

############################ construction of the algorithm

# helpers
powerMat <- function(m, power = 1) {
  if (power == 0) {
    return(diag(1, ncol(m)))
  } else {
    return(m%*%powerMat(m, power = power-1))
  }
}


# chain
p.random.nei <- 0.05 # probability of moving to a random neighbor node
p.staying <- 0.05 # probability of staying at the current node

# randomize variables and splits
variable.mat <- graph.mat
for (i in 1:nrow(graph.mat)) {
  variable.mat[i, graph.mat[i,]==1] <- sample(1:d, 
                                             1, 
                                       replace=TRUE)
}
split.mat <- graph.mat
for (i in 1:nrow(graph.mat)) {
  v <- variable.mat[i,graph.mat[i,]==1][1]
  split.mat[i, graph.mat[i,]==1] <- sample(X[,v],1)
}
split.mat[graph.mat==1] <- sapply(variable.mat[graph.mat==1], function(v) sample(X[,v],1))



piStationary <- function(x, method = "eigen") {
  
  # if we have missing values, we sub with Inf
  x[!is.finite(x)] <- NA
  
  p <- sapply(1:sum(graph.mat==1), 
              function(j) 
                as.numeric(x[variable.mat[graph.mat==1][j]]<=split.mat[graph.mat==1][j]))
  
  pM <- graph.mat
  for (i in 1:nrow(graph.mat)) {
    xval <- x[variable.mat[i, graph.mat[i,]==1][1]]
    cond <- xval <= split.mat[i,graph.mat[i,]==1][1]
    pM[i, graph.mat[i,]==1][1] <- (1-p.random.nei)*ifelse(is.na(xval), 1/2, as.numeric(cond)) + p.random.nei * 1/2
    pM[i, graph.mat[i,]==1][2] <- (1-p.random.nei)*ifelse(is.na(xval), 1/2, 1-as.numeric(cond)) + p.random.nei * 1/2
  }
  
  prob <- (p.staying)*diag(1,nrow(pM),ncol(pM)) + (1-p.staying)*pM
  #prob <- pM
  # raise to the power
  probP <- powerMat(prob,1000)
  
  if (method == "eigen") {
    eigen(t(prob))$vectors[,1]
    pistat <- abs(eigen(t(prob))$vectors[,1]) / sum(abs(eigen(t(prob))$vectors[,1]))
  } else {
    probP <- powerMat(t(prob),300)
    pistat <- apply(probP,1,mean)
  }
  
  return(pistat)
}




############################ testing
piMat_na <- t(apply(Xna, 1, piStationary))
piMat <- t(apply(X[,1:4], 1, piStationary))

# testing
Xtest_na <- Xtest
Xtest_na[,5] <- NA

piMat_test <- t(apply(Xtest, 1, piStationary))
piMat_test_na <- t(apply(Xtest_na, 1, piStationary))

k <- piMat %*% t(piMat)
kna <- piMatna %*% t(piMatna)


# setup no NA at all!!!!

# linear regression prediction from X1-X4 on X5
l <- lm(x.5~., data = data.frame(x=X))
preds_lm <- predict(l, data.frame(x=Xtest))

# ranger prediction from X1-X4 on X5
rf <- ranger::ranger(x.5~., data = data.frame(x=X), mtry = 4, num.trees = 2000)
preds_rf <- predict(rf, data.frame(x=Xtest))$predictions

# our approach
lpi <- lm(y~., data = data.frame(y=X[,5],p=piMat))
preds_lm_pi <- predict(lpi, data.frame(p=piMat_test))
preds_lm_pi_na <- predict(lpi, data.frame(p=piMat_test_na))

# MSE's
mean((mean(Xtest[,5]) - Xtest[,5])^2)
mean((preds_lm - Xtest[,5])^2)
mean((preds_rf - Xtest[,5])^2)
mean((preds_lm_pi - Xtest[,5])^2)
mean((preds_lm_pi_na - Xtest[,5])^2)


# - ensemble
# - data driven selection of splits?
# - range of applications?
# - tuning parameters (mixture of possibilities)
# - assumptions for the method to work