## Script of Nicolai, 2.4.20


library(randomForest)

## generate inverse-smiley data
n=1000
x = runif(n,-1,1)
X = cbind(x, -x^2+1 + rnorm(n)*0.02)
plot(X)

## sample from a null distribution (destroying the copula)
Z = cbind(X[,1], sample(X[,2], n))

## classify between true and null distribution to get density estimate with RF
XZ = rbind(X,Z)
colnames(X) = colnames(Z) = colnames(XZ) = paste("Var", 1:ncol(X),sep="")
Y = c(rep(1, nrow(X)), rep(0,nrow(Z)))

rf = randomForest(XZ, as.factor(Y),keepnodes=TRUE)
print(rf)

## check which nodes of the RF contain true and null samples
nodesX = attr(predict( rf, X, nodes=TRUE),"nodes")
#nodesZ = attr(predict( rf, Z, nodes=TRUE), "nodes")

## predict missing value of X_1 if X_2 = 0.5
x2 = 0.5
tryat = cbind( sample(X[,1],n), rep(x2, n))
colnames(tryat) = colnames(X)
nodestry = attr( predict(rf, tryat, nodes=TRUE),"nodes")
inXnodes  = 0*nodestry
for (k in 1:ncol(nodestry)){
  inXnodes[,k]  = as.numeric( nodestry[,k] %in% nodesX[,k])
}
weights = apply(inXnodes,1,sum)
plot(tryat[,1], weights)
## where to cut here to get CI?? probably somewher close to 450 out of
# 500 trees or thereabouts...

############### simple rf ######################################################
#### imbalanced? more weights on the positive side...
dat <- data.frame(x1 = X[,1], x2 = X[,2])

rf <- ranger::ranger(formula = x1~x2,data=dat, quantreg=T, num.trees = 5000, 
                     splitrule = "extratrees")
imp <- predict(rf, data= data.frame(x2 = 0.5), type="quantiles",
                     quantiles = seq(0.1,0.9, length.out = 1000))$predictions

hist(c(imp),breaks=10)
length(which(c(imp)<=-0.5))/length(c(imp))
length(which(c(imp)>=0.5))/length(c(imp))

mean(dat$x1>=0)


#### QUESTIONS:
## what is the probability that a sample from the true distribution
## would fall into the observed nodes nodesX in at least 1-alpha of all
## trees (determines cutoff needed for CI) ?
  
## for more noisy data: might need to exclude nodes in nodesX in which
## only very few true observations fall and accordingly lower the
## thresshold for the number of trees necessary

## if more than one variable is missing it might be more interesting.
## if we have many trees, might be able to get CI for *some* of the
## missing variables; can then enter these bounds for the missing
## variables and check the possible nodes again, getting tighter bounds
## and repeat until convergence. needs careful analysis of error
## probabilities to make sure we get valid CIs in the end.