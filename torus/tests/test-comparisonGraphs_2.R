

# use torus package
require(torus)
n<-500
d <- genData(n =n , dataset = "bivariateGaussian", pNA = 0)

# plot the data
plot(d$X,pch=19)

# grid of interesting points
X.NA.grid <- matrix(c(NA, -20,
                      NA, -12,
                      NA, -8,
                      NA, -4,
                      NA, 0),ncol=2,byrow = T)



ct<-c()
et<-c()
B<-1

if (B > 1) {

for (b in 1:B){

  ct[[b]]<-closedTree(X = d$X.NA, depth = 3)
  et[[b]] <- extraTorus(X = d$X.NA, nb.nodes = 10)

}


w_ct_all <- getSampleWeightsEnsemble(ct, X = rbind(X.NA.grid,d$X.NA),
                             prob = c(0.005,0.005,1-2*0.005), d="k_divergence", method="eigen")
w_et <- getSampleWeightsEnsemble(et, X = rbind(X.NA.grid,d$X.NA),
                                     prob = c(0.005,0.005,1-2*0.005), d="k_divergence", method="eigen")
}else{



et <- extraTorus(X = d$X.NA, nb.nodes = 10)
w_et<-getSampleWeights(et, X = rbind(X.NA.grid,d$X.NA),
                                        prob = c(0.005,0.005,1-2*0.005), d="k_divergence", method="eigen" ) # ,  method = "power", power=2

# around 2 min for both methods





ct <- closedTree(X = d$X.NA, depth = 4)
w_ct_all<-getSampleWeights(ct, X = rbind(X.NA.grid,d$X.NA),
                       prob = c(0.005,0.005,1-2*0.005), d="k_divergence", method="eigen")

}


if (  sum(is.na(d$X.NA))==0 ){


  ## If we have no missing values, check nearest neigbours of the first point
par(mfrow=c(2,1))
plot(d$X.NA)
points(d$X.NA[order(w_et[6,7:(n+5)], decreasing = FALSE)[1:round(n/5)],], col="green")
points(d$X.NA[1,1],d$X.NA[1,2], col="red" )


plot(d$X.NA)
points(d$X.NA[order(w_ct_all[6,7:(n+5)], decreasing = FALSE)[1:round(n/5)],], col="green")
points(d$X.NA[1,1],d$X.NA[1,2], col="red" )
}

par(mfrow=c(1,1))
plot(d$X.NA, cex=w_et[6,7:(n+5)]*100, pch=19)
points(d$X.NA[1,1],d$X.NA[1,2], col="red" )





par(mfrow=c(5,2))
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19, main=toString(X.NA.grid[i,]))
  points(d$X.NA[order(w_et[i,6:(n+5)], decreasing = FALSE)[1:round(n/10)],],col="green",pch=19)
  points(d$X.NA[order(w_et[i,6:(n+5)], decreasing = FALSE)[1:round(n/50)],],col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),pch=19)
  points( matrix(c(mean(w_et[i,6:(n+5)]*d$X.NA[,1], na.rm = T), X.NA.grid[i,2]), ncol = 2)  ,col="red",pch=19)

  plot(d$X.NA,pch=19, main=toString(X.NA.grid[i,]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = FALSE)[1:round(n/10)],],col="green",pch=19)
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = FALSE)[1:round(n/100)],],col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),pch=19)
  points( matrix(c(mean(w_ct_all[i,6:(n+5)]*d$X.NA[,1], na.rm = T), X.NA.grid[i,2]), ncol = 2)  ,col="red",pch=19)
}


# If there are no NA we can compare with
distmat<-distance(d$X.NA , method="euclidean")












