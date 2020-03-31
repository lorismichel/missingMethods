genData <- function(n, dataset = "bivariateGaussian", pNA = 0.05) {
  
  if (dataset == "bivariateGaussian") {
    # data
    n <- 500 # number of observations
    d <- 2 # number of dimensions
    
    X <- MASS::mvrnorm(n = n, mu = rep(0, d),
                       Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example
    Xtest <- MASS::mvrnorm(n = n, mu = rep(0, d),
                           Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example
    
    # NA filter
    na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(1-pNA, pNA), replace = TRUE))
    X.NA <- X
    X.NA[na==1] <- NA
    
    return(list(X = X, X.NA = X.NA))
   
  } else {
    
    d <- 2 
    
    X1 <- rnorm(n)
    X2<- -2*X1^2 + rnorm(n)*0.5
    
    X <- cbind(X1,X2)
   
    
    
    # NA filter
    na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(1-pNA, pNA), replace = TRUE))
    X.NA <- X
    X.NA[na==1] <- NA
    
    return(list(X = X, X.NA = X.NA))
  }
  
  
  
}

