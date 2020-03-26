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

