

# use torus package
require(torus)
n<-1000
d <- genData(n =n , dataset = "parabola", pNA = 0.05)

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
B<-20
for (b in 1:B){

  ct[[b]]<-closedTree(X = d$X.NA, depth = 4)
  et[[b]] <- extraTorus(X = d$X.NA, nb.nodes = 4)

}


w_ct_all <- getSampleWeightsEnsemble(ct, X = rbind(X.NA.grid,d$X.NA),
                             prob = c(0.005,0.005,1-2*0.005), method="eigen")
w_et <- getSampleWeightsEnsemble(et, X = rbind(X.NA.grid,d$X.NA),
                                     prob = c(0.005,0.005,1-2*0.005), method="eigen")



par(mfrow=c(5,2))
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19, main=toString(X.NA.grid[i,]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_et[i,], decreasing = TRUE)[1:round(n/10)],],col="green",pch=19)
  points(rbind(X.NA.grid,d$X.NA)[order(w_et[i,], decreasing = TRUE)[1:round(n/100)],],col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),pch=19)
  points( matrix(c(mean(w_et[i,6:(n+5)]*d$X.NA[,1], na.rm = T), X.NA.grid[i,2]), ncol = 2)  ,col="red",pch=19)

  plot(d$X.NA,pch=19, main=toString(X.NA.grid[i,]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = TRUE)[1:round(n/10)],],col="green",pch=19)
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,], decreasing = TRUE)[1:round(n/100)],],col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),pch=19)
  points( matrix(c(mean(w_ct_all[i,6:(n+5)]*d$X.NA[,1], na.rm = T), X.NA.grid[i,2]), ncol = 2)  ,col="red",pch=19)
}


par(mfrow=c(5,2))
## looking at the points
for (i in 1:nrow(X.NA.grid)) {
  plot(d$X.NA,pch=19, main=toString(X.NA.grid[i,]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_et[i,w_et[i,] - 1/n > 0], decreasing = T),],col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),pch=19)
  points( matrix(c(mean(w_et[i,6:(n+5)]*d$X.NA[,1], na.rm = T), X.NA.grid[i,2]), ncol = 2)  ,col="red",pch=19)

  plot(d$X.NA,pch=19, main=toString(X.NA.grid[i,]))
  points(rbind(X.NA.grid,d$X.NA)[order(w_ct_all[i,w_ct_all[i,] - 1/n > 0], decreasing = T),],col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),pch=19)
  points( matrix(c(mean(w_ct_all[i,6:(n+5)]*d$X.NA[,1], na.rm = T), X.NA.grid[i,2]), ncol = 2)  ,col="red",pch=19)
}








