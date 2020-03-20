# Random node graph 

# architecture of graph nodes
#1 -> 2 -> 3
#-    |    -
#|    -    |
#4 <- 5 -> 6
#|    -    |
#-    |    -
#7 -> 8 <- 9

# define global parameters

# data
n <- 200 # number of observations
d <- 5 # number of dimensions

X <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = diag(0.6,d,d)+matrix(0.4,d,d)) # naive matrix example
Xtest <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = diag(0.6,d,d)+matrix(0.4,d,d)) # naive matrix example

# NA filter
na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(0.95, 0.05), replace = TRUE))
Xna <- X
Xna[na==1] <- NA
# chain
m <- 0.01 # probability of moving to a random neighbor node

# graph connexions
graph.mat <- matrix(0,ncol=9,nrow=9)
graph.mat[5,4] <- 1
graph.mat[4,1] <- 1
graph.mat[1,2] <- 1
graph.mat[2,5] <- 1
graph.mat[5,6] <- 1
graph.mat[6,3] <- 1
graph.mat[3,2] <- 1
graph.mat[4,7] <- 1
graph.mat[7,8] <- 1
graph.mat[8,5] <- 1
graph.mat[6,9] <- 1
graph.mat[9,8] <- 1

# power mat
powerMat <- function(m, power = 1) {
  if (power == 0) {
    return(diag(1, ncol(m)))
  } else {
    return(m%*%powerMat(m, power = power-1))
  }
}

# randomize variables and splits
variable.mat <- graph.mat
variable.mat[graph.mat==1] <- sample(1:d, sum(graph.mat==1), replace=TRUE)

split.mat <- graph.mat
split.mat[graph.mat==1] <- sapply(variable.mat[graph.mat==1], function(v) sample(X[,v],1))

piStationary <- function(x) {
  x[!is.finite(x)] <- Inf

  p <- sapply(1:sum(graph.mat==1), function(j) as.numeric(x[variable.mat[graph.mat==1][j]]<=split.mat[graph.mat==1][j]))
  pM <- graph.mat
  pM[graph.mat==1] <- p
  pM <- t(apply(pM,1,function(x) x/max(1,sum(x))))
  di <- apply(pM,1,function(x) 1-sum(x))
  diag(pM) <- diag(pM) + di
  
  pMc <- graph.mat
  pMc <- m*pMc
  #pMc <- t(apply(pMc,1,function(x) x/max(1,sum(x))))
  diag(pMc) <- diag(pMc) + apply(pMc,1,function(x) 1-sum(x))
  
  prob <- (pM + pMc)/2
  probP <- powerMat(prob,300)
  return(apply(probP,2,mean))
}
piMatna <- t(apply(Xna, 1, piStationary))
piMat <- t(apply(X, 1, piStationary))
piMat_test <- t(apply(Xtest, 1, piStationary))

k <- piMat %*% t(piMat)
kna <- piMatna %*% t(piMatna)


# prediction of last column
l <- lm(x.5~., data = data.frame(x=X))
preds_lm <- predict(l, data.frame(x=Xtest))

lpi <- lm(y~., data = data.frame(y=X[,5],p=piMat))
preds_lm_pi <- predict(lpi, data.frame(p=piMat_test))

mean((preds_lm-Xtest[,5])^2)
mean((preds_lm_pi-Xtest[,5])^2)
mean((mean(Xtest[,5])-Xtest[,5])^2)


