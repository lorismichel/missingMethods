############################ testing extraTorus

# use torus package
require(torus)
# data
n <- 500 # number of observations
d <- 2 # number of dimensions

X <- MASS::mvrnorm(n = n, mu = rep(0, d),
                   Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example
Xtest <- MASS::mvrnorm(n = n, mu = rep(0, d),
                       Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example

# NA filter
na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(0.95, 0.05), replace = TRUE))
X.NA <- X
X.NA[na==1] <- NA
X.NA.grid <- matrix(c(-2, NA, 2, NA, NA, 0),ncol=2,byrow = T)

# plot the data
plot(X,pch=19)

# generate an extraTorus
tor1 <- extraTorus(X = X.NA, nb.nodes = 4)
tor2 <- extraTorus(X = X.NA, nb.nodes = 4)
tor3 <- extraTorus(X = X.NA, nb.nodes = 4)
tor4 <- extraTorus(X = X.NA, nb.nodes = 4)

to <- combine(tor1, tor2)
# get weights from stationary distr
weights <- getSampleWeights(to, X = rbind(X.NA.grid,X.NA))

# vizu results
par(mfrow=c(1,2))
plot(X,pch=19)
plot(X.NA,pch=19)

# investigate first two points
plot(rbind(X.NA.grid,X.NA), cex = weights[1,]*500,pch=19)
plot(rbind(X.NA.grid,X.NA), cex = weights[2,]*500,pch=19)

spatstat::weighted.quantile(x = rbind(X.NA.grid,X.NA)[,2], w = weights[1,], probs = c(0.05,0.95))
spatstat::weighted.quantile(x = rbind(X.NA.grid,X.NA)[,2], w = weights[2,], probs = c(0.05,0.95))





######## nonlinear examples


n <- 1000 # number of observations
d <- 2 # number of dimensions

X1 <- rnorm(n)
X2<- -2*X1^2 + rnorm(n)*0.5

X=cbind(X1,X2)



# NA filter
na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(0.95, 0.05), replace = TRUE))
X.NA <- X
X.NA[na==1] <- NA
X.NA.grid <- matrix(c(-2, NA, 2, NA, NA, 0),ncol=2,byrow = T)


# generate an extraTorus
to<-c()
B<-200
for (b in 1:B){

  to[[b]]<-extraTorus(X = X.NA, nb.nodes = 5)

}


# get weights from stationary distr
weights <- getSampleWeightsEnsemble(to, X = rbind(X.NA), method="eigen")

# vizu results
X.NAtrue<-X[!complete.cases(X.NA),]
dim(X.NAtrue)[2]

par(mfrow=c(1,3))
plot(X.NA ,pch=19)
points(X.NAtrue,pch=19, col="blue")

# investigate imputed points (Jeff)
X.indexed<-cbind(1:n, X.NA)
X.imputed<- apply(X.indexed[!complete.cases(X.indexed),],1, function(x){ colSums(weights[x[1],]*X.NA, na.rm = T) } )
X.imputed <- t(X.imputed)
plot(X.NA ,pch=19)
points(X.imputed,pch=19, col="red")

X.imputed <- X.NA
inds_missing <-c()
# imputation by ML
for (j in 1:d){
  inds_to_impute <- which(is.na(c(X.NA[,j]))==TRUE)
  for(k in inds_to_impute){
    X.imputed[k,j] <- colSums(weights[k,]*X.NA, na.rm = T)[j] #mean(weights[k,]*X.NA[,j], na.rm=T)
  }
  inds_missing <- c(inds_to_impute, inds_missing)
}

plot(X.NA, pch=19)
points(X.imputed[inds_missing,], col="blue", pch=19)





## Investigate distribution for a single torus
pi1=stationaryDistr(object = tor1, X = X.NA, method="eigen")
par(mfrow=c(3,3))


for (i in 1:9) {

  barplot(pi1[i,], ylim=c(0,0.5), main=paste( toString(round(X[i,],3)))  )


}


if (type == "nested") {
  
  # compute the parent node vector (only used for non-terminal node)
  parent.nodes.vec <- apply(adj.mat, 2, function(x) which(x==1)[1])
  
  # this function works within a tree
  getSamplesNode <- function(node, X) {
    
    # the direct parent node
    parent.node <- parent.nodes.vec[node]
    
    # stopping condition
    if (node == 1) {
      
      return(1:nrow(X))
      
    } else {
      
      sub <- which(apply(X, 1, function(x) {
        
        # does it go left or right?
        dec <- which(node == which(ct1$adj.mat[parent.node,]==1))
        
        # decide
        if (dec == 1) {
          x[variable.mat[parent.node, which(adj.mat[parent.node,]==1)[1]]] <= split.mat[parent.node, which(adj.mat[parent.node,]==1)[1]]
        } else {
          x[variable.mat[parent.node, which(adj.mat[parent.node,]==1)[1]]] > split.mat[parent.node, which(ct1$adj.mat[parent.node,]==1)[1]]
        }
      }))
    }
    
    return(intersect(sub, getSamplesNode(parent.node, X)))
  } 
  
  
  
  
  
  
}




