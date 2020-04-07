genData <- function(n, d=2, dataset = "bivariateGaussian", pNA = 0.05) {

  if (dataset == "bivariateGaussian") {
    X <- MASS::mvrnorm(n = n, mu = rep(0, d),
                       Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example

    # NA filter
    na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(1-pNA, pNA), replace = TRUE))
    X.NA <- X
    X.NA[na==1] <- NA

    return(list(X = X, X.NA = X.NA))

  } else {

    x = runif(n,-1,1)
    X = cbind(x, -x^2+1 + rnorm(n)*0.02)

    # NA filter
    na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(1-pNA, pNA), replace = TRUE))
    X.NA <- X
    X.NA[na==1] <- NA

    return(list(X = X, X.NA = X.NA))
  }
}

