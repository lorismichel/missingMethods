############################ testing extraTorus

# data
n <- 500 # number of observations
d <- 2 # number of dimensions

X <- MASS::mvrnorm(n = n, mu = rep(0, d),
                   Sigma = diag(0.6,d,d)+matrix(0.4,d,d)) # naive matrix example
Xtest <- MASS::mvrnorm(n = n, mu = rep(0, d),
                       Sigma = diag(0.6,d,d)+matrix(0.4,d,d)) # naive matrix example

# NA filter
na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(0.95, 0.05), replace = TRUE))
X.NA <- X
X.NA[na==1] <- NA


# generate an extraTorus
tor <- extraTorus(X = X, nb.nodes = 2)

# get the stationary distribution
pi.power <- stationaryDistr(object = tor, X.NA, method = "power")
pi.eigen <- stationaryDistr(object = tor, X.NA, method = "eigen")


# get weights for first observation
k <- apply(pi.eigen, 1, function(x) sum(pi.eigen[500,]*x))
k <- k/sum(k)

plot(X,pch=19)
plot(X, cex = k*500,pch=19)
points(X[500,1],X[500,2],col="red",pch=19,cex=3)

X.NA[500,]
X[500,]

pX <- X.NA[sample(1:nrow(X), prob = k, size = 10000, replace = TRUE),]
mean(pX[,1],na.rm = TRUE)
quantile(pX[,1],na.rm = TRUE, probs = 0.025)
quantile(pX[,1],na.rm = TRUE, probs = 0.975)



plot(X,pch=19)
points(pX,pch=19,col="blue")


plot(density(na.omit(pX[,1])))
# 1) very naive baseline:
# for each observation that we have, sample according to stationary a random node --> forest game
# for a given point x --> pi.stat --> resample a node based on pi.stat
# down-side: for large graphs it is very variance prone --> correction allowing multiple resampling

# 2) kernel analogy with random forests leaf node
# if we have a tree, the empirical distribution is x=(1,0,0)
# (1,0,0) * (0,1,0), (1,0,0) * (1,0,0) only one leaf node (point mass)
# w_i(x) = pi_x * pi_x_i / sum(p_x * p_x_i) smooth assignement

# brainstorm on how to impute:
# - i take the folowing mechanism: for point x_1 and x_2 select according to pi_1 and pi_2 the node independly
# - what is the probability of both of them falling into the same node. c(0.1,0.2,0.7) c(0,0.5, 0.5) propto 0.2 * 0.5 + 0.5 * 0.7 =


