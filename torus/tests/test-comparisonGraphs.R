############################ comparing different graphs

####################################################################################################
################################ example 1 #########################################################
####################################################################################################
# Simple example, using just one torus/tree, no ensemble
# use torus package
require(torus)

d <- genData(n = 800, d=2, dataset = "parabola", pNA = 0.05)

# plot the data
plot(d$X,pch=19, xlab="X1", ylab="X2")

# grid of interesting points

X.NA.grid <- matrix(c(NA, 0,
                      NA, 0.2,
                      NA, 0.5,
                      NA, 0.8,
                      NA, 1),ncol=2,byrow = T)



# generate an extraTorus and a closed tree
et <- extraTorus(X = rbind(X.NA.grid,d$X.NA), nb.nodes = 10)
ct <- closedTree(X = rbind(X.NA.grid,d$X.NA), depth = 4)

# look at transition matrix: how close are they?
tr1 <- getTransitionMatrix(ct, x = c(NA, 0))
tr2 <- getTransitionMatrix(ct, x = c(NA, 0.2))
mean(abs(tr1-tr2))


# get weights from stationary distr

prob <- c(0.0001,0.0001,1-2*0.0001)
st_et <- stationaryDistr(et, rbind(X.NA.grid,d$X.NA), prob = prob, method = "eigen")
st_ct_all <- stationaryDistr(ct, rbind(X.NA.grid,d$X.NA), prob = prob, method = "eigen")
st_ct_leaves  <- stationaryDistr(ct, rbind(X.NA.grid,d$X.NA), prob = prob, subset = TRUE, method = "eigen")


############### looking at stationary distributions ###############################################
####### stationary distr for one torus/tree

par(mfrow=c(5,3))
for (i in 1:nrow(X.NA.grid)) {
  plot(st_et[i,],pch=19, main =paste0("Extra Torus, X2 = ", X.NA.grid[i,2]))
  plot(st_ct_all[i,],pch=19, main =paste0("Closed Tree ALl, X2 = ", X.NA.grid[i,2]))
  plot(st_ct_leaves[i,],pch=19, main = paste0("Closed Tree Leafs, X2 = ", X.NA.grid[i,2]))
}


############### looking at sample weights #########################################################

prob <- c(0.0001,0.0001,1-2*0.0001)
w_et <- getSampleWeights(et, X = rbind(X.NA.grid,d$X.NA),
                            prob = prob, method = "eigen", power=100)

w_ct_all <- getSampleWeights(ct, X = rbind(X.NA.grid,d$X.NA),
                             prob = prob, method = "eigen")
w_ct_leaves <- getSampleWeights(ct, X = rbind(X.NA.grid,d$X.NA),
                                prob = prob,
                            subset = TRUE, method = "eigen")

par(mfrow=c(5,3))
nr_neighbours <- 100
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19, main =paste0("Extra Torus, X2 = ", X.NA.grid[i,2]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_et[i,], decreasing = TRUE)[1:nr_neighbours],],col="red",pch=19)
  plot(d$X.NA,pch=19, main = paste0("Closed Tree All, X2 = ", X.NA.grid[i,2]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = TRUE)[1:nr_neighbours],],col="red",pch=19)
  plot(d$X.NA,pch=19, main = paste0("Closed Tree Leafs, X2 = ", X.NA.grid[i,2]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:nr_neighbours],],col="red",pch=19)
}

######################################################################################################
############################## Forest Approach #######################################################
######################################################################################################
### Can we see the modes in easy random forest usage?

# forest approach
rf <- ranger::ranger(formula = X1~X2, data = data.frame(X1=na.omit(d$X.NA)[,1],X2=na.omit(d$X.NA)[,2]),
                     splitrule =  "extratrees", quantreg = TRUE, num.random.splits = 1)
preds <- predict(rf, data.frame(X2=X.NA.grid[,2]), type = "quantiles", quantiles = runif(1000))$predictions


par(mfrow=c(5,2))
for (i in 1:nrow(X.NA.grid)) {
  plot(density(na.omit(rbind(X.NA.grid,d$X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:100],1])),
       main=paste0("Weights of Leaves, X2 = ",X.NA.grid[i,2]),xlab="X1",ylab="X2",xlim=c(-4,4))
  plot(density(preds[i,]),xlab="X1",ylab="X2",main=paste0("Quantiles of RF, X2 = ",X.NA.grid[i,2]),xlim=c(-4,4))
}

######################################################################################################
############################### Example 2 ############################################################
######################################################################################################
# Ensemble
d <- genData(n = 300, d=2, dataset = "parabola", pNA = 0.05)
et <- NULL
ct <- NULL
B <- 5
for (b in 1:B){
  ct[[b]]<- closedTree(X = d$X.NA, depth = 5)
  et[[b]] <- extraTorus(X = d$X.NA, nb.nodes =12)
}

w_et <- getSampleWeightsEnsemble(et, X = rbind(X.NA.grid,d$X.NA),
                         prob = c(0.05,0.05,1-2*0.05))$w
w_ct_all <- getSampleWeightsEnsemble(ct, X = rbind(X.NA.grid,d$X.NA),
                                 prob = c(0.05,0.05,1-2*0.05), subset = FALSE)$w
w_ct_leaves <- getSampleWeightsEnsemble(ct, X = rbind(X.NA.grid,d$X.NA),
                                 prob = c(0.05,0.05,1-2*0.05), subset = TRUE)$w


par(mfrow=c(5,3))
## looking at the sample weight points
nr_neighbours <- 80
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19,main = paste0("Extra Torus, X2 = ", X.NA.grid[i,2]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_et[i,], decreasing = TRUE)[1:nr_neighbours],],col="red")
  plot(d$X.NA,pch=19,main = paste0("Closed Tree All, X2 = ", X.NA.grid[i,2]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = TRUE)[1:nr_neighbours],],col="red")
  plot(d$X.NA,pch=19,main = paste0("Closed Tree Leaves, X2 = ", X.NA.grid[i,2]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_leaves[i,], decreasing = TRUE)[1:nr_neighbours],],col="red")
}

#### looking at the stationary distr
pi_et <- getSampleWeightsEnsemble(et, X = rbind(X.NA.grid,d$X.NA),
                                 prob = c(0.05,0.05,1-2*0.05))$pi.mat
pi_ct_all <- getSampleWeightsEnsemble(ct, X = rbind(X.NA.grid,d$X.NA),
                                     prob = c(0.05,0.05,1-2*0.05), subset = FALSE)$pi.mat
pi_ct_leaves <- getSampleWeightsEnsemble(ct, X = rbind(X.NA.grid,d$X.NA),
                                        prob = c(0.05,0.05,1-2*0.05), subset = TRUE)$pi.mat

par(mfrow=c(5,3))
for (i in 1:nrow(X.NA.grid)) {
  plot(pi_et[i,],pch=19,main = paste0("Extra Torus, X2 = ", X.NA.grid[i,2]))
  plot(pi_ct_all[i,],pch=19,main = paste0("Closed Tree All, X2 = ", X.NA.grid[i,2]))
  plot(pi_ct_leaves[i,],pch=19,main = paste0("Closed Tree Leaves, X2 = ", X.NA.grid[i,2]))
}



##################################### what are these functions?##################################
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


################################ What is tried here? ###############################################

getSamplesNode(node = 6, X = d$X)

get


tt <- ranger::ranger(formula = Y~X1+X2, data = data.frame(Y=rnorm(nrow(na.omit(d$X.NA))),
                                                          X1=na.omit(d$X.NA)[,1],X2=na.omit(d$X.NA)[,2]),
                     splitrule =  "extratrees", quantreg = TRUE, num.random.splits = 1, max.depth = 3, num.trees = 1)

tt$forest$split.varIDs
tt$forest$split.values
tt$forest$child.nodeIDs
