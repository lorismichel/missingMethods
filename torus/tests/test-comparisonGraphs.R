############################ comparing different graphs

# use torus package
require(torus)

d <- genData(n = 1000, dataset = "parabola", pNA = 0.05)

# plot the data
plot(d$X,pch=19)

# grid of interesting points
X.NA.grid <- matrix(c(NA, -20, 
                      NA, -12, 
                      NA, -8, 
                      NA, -4, 
                      NA, -1,
                      NA, 0),ncol=2,byrow = T)



# generate an extraTorus and a closed tree
et <- extraTorus(X = d$X.NA, nb.nodes = 40)
ct <- closedTree(X = d$X.NA, depth = 4)

# look at transition matrix
tr1 <- getTransitionMatrix(ct, x = c(NA,-20))
tr2 <- getTransitionMatrix(ct, x = c(NA,-4))
mean(abs(tr1-tr2))



# get weights from stationary distr
id.leaves <- which(apply(ct$adj.mat,1,sum)==1)

prob <- c(0.0001,0.0001,1-2*0.0001)
st_et <- stationaryDistr(et, rbind(X.NA.grid,d$X.NA), prob = prob, method = "eigen")
st_ct_all <- stationaryDistr(ct, rbind(X.NA.grid,d$X.NA), prob = prob, method = "eigen")
st_ct_leaves  <- stationaryDistr(ct, rbind(X.NA.grid,d$X.NA), prob = prob, subset = id.leaves, method = "eigen")


par(mfrow=c(6,3))
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(st_et[i,],pch=19)
  plot(st_ct_all[i,],pch=19)
  plot(st_ct_leaves[i,],pch=19)
}


getTransitionMatrix(object = ct, x = c(NA,-2))
getTransitionMatrix(object = ct, x = c(NA,2))


par(mfrow=c(3,1))

ct <- closedTree(X = d$X.NA, depth = 4)
et <- extraTorus(X = d$X.NA, nb.nodes = 4)

prob <- c(0.0001,0.0001,1-2*0.0001)
w_et <- getSampleWeights(et, X = rbind(X.NA.grid,d$X.NA), 
                            prob = prob, method = "power", power=100)
w_ct_all <- getSampleWeights(ct, X = rbind(X.NA.grid,d$X.NA),
                             prob = prob, method = "eigen")
w_ct_leaves <- getSampleWeights(ct, X = rbind(X.NA.grid,d$X.NA), 
                                prob = prob,
                            subset = id.leaves, method = "eigen")

par(mfrow=c(6,3))
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19)
  points(rbind(X.NA.grid,d$X.NA)[order(w_et[i,], decreasing = TRUE)[1:100],],col="red",pch=19)
  plot(d$X.NA,pch=19)
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = TRUE)[1:100],],col="red",pch=19)
  plot(d$X.NA,pch=19)
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:100],],col="red",pch=19)
}



# forest approach
rf <- ranger::ranger(formula = X1~X2, data = data.frame(X1=na.omit(d$X.NA)[,1],X2=na.omit(d$X.NA)[,2]),
                     splitrule =  "extratrees", quantreg = TRUE, num.random.splits = 1)
preds <- predict(rf, data.frame(X2=X.NA.grid[,2]), type = "quantiles", quantiles = runif(1000))$predictions

par(mfrow=c(2,3))


rbind(X.NA.grid,d$X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:100],2]

par(mfrow=c(6,2))
for (i in 1:nrow(X.NA.grid)) {
  plot(density(na.omit(rbind(X.NA.grid,d$X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:100],1])),main="",xlab="X1",ylab="X2",xlim=c(-4,4))
  plot(density(preds[i,]),xlab="X1",ylab="X2",main="",xlim=c(-4,4))
}



ct1 <- closedTree(X = d$X.NA, depth = 2)
ct2 <- closedTree(X = X.NA, depth = 3)
ct3 <- closedTree(X = X.NA, depth = 3)
ct4 <- closedTree(X = X.NA, depth = 3)
ct5 <- closedTree(X = X.NA, depth = 3)

f <- combine(ct1,ct2,ct3,ct4,ct5)
w_f <- getSampleWeightsEnsemble(f, X = rbind(X.NA.grid,X.NA), 
                         prob = c(0.005,0.005,1-2*0.005), subset = id.leaves)



par(mfrow=c(5,3))
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19)
  points(rbind(X.NA.grid,X.NA)[order(w_f[i,], decreasing = TRUE)[1:200],],col="red")
  plot(d$X.NA,pch=19)
  points(rbind(X.NA.grid,X.NA)[order(w_ct_all[i,], decreasing = TRUE)[1:200],],col="red")
  plot(d$X.NA,pch=19)
  points(rbind(X.NA.grid,X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:200],],col="red")
}


dropDown <- function(ct, x) {
  path <- c()
  nb.children <- apply(ct$adj.mat,1,sum)
  cur.node <- 1
  
  while(nb.children[cur.node]!=1) {
    path <- c(path,cur.node)
    next.nodes.candidates <- which(ct$adj.mat[cur.node,]==1)
    x[ct$variable.mat[i,next.nodes.candidates][1]] <= ct$split.mat[i,next.nodes.candidates[1]]
  }
  
  
}

parent.nodes.vec <- apply(ct1$adj.mat, 2, function(x) which(x==1)[1])

# this function works within a tree
getSamplesNode <- function(node, X) {
  
  # the direct parent node
  parent.node <- parent.nodes.vec[node]
  #parent.parent.node <- parent.nodes.vec[parent.node]
  if (node == 1) {
    return(1:nrow(X))
  } else {
    sub <- which(apply(X, 1, function(x) {
      dec <- which(node == which(ct1$adj.mat[parent.node,]==1))
      if (dec == 1) {
        x[ct1$variable.mat[parent.node,which(ct1$adj.mat[parent.node,]==1)[1]]] <= ct1$split.mat[parent.node, which(ct1$adj.mat[parent.node,]==1)[1]]
      } else {
        x[ct1$variable.mat[parent.node,which(ct1$adj.mat[parent.node,]==1)[1]]] > ct1$split.mat[parent.node, which(ct1$adj.mat[parent.node,]==1)[1]]
      }
    }))
  }
  
  return(intersect(sub, getSamplesNode(parent.node, X)))
} 




getSamplesNode(node = 6, X = d$X)

get


tt <- ranger::ranger(formula = Y~X1+X2, data = data.frame(Y=rnorm(nrow(na.omit(d$X.NA))), X1=na.omit(d$X.NA)[,1],X2=na.omit(d$X.NA)[,2]),
                     splitrule =  "extratrees", quantreg = TRUE, num.random.splits = 1, max.depth = 3, num.trees = 1)

tt$forest$split.varIDs
tt$forest$split.values
tt$forest$child.nodeIDs
