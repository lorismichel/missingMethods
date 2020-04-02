powerMat <- function(m, power = 1) {
  if (power == 0) {
    return(diag(1, ncol(m)))
  } else {
    return(expm::`%^%`(m, power))
  }
}
