
SimuData <- function(Loadings, FactorCov, N){
  CovMatrix <- Loadings %*% FactorCov %*% t(Loadings)
  diag(CovMatrix) <- 1
  if (any(eigen(CovMatrix)$values < 0)){
    return(print("The population variance covariance matrix are not positive defined"))
  }
  Means <- rep(0, each = ncol(CovMatrix))
  simu_data <-data.frame(mvrnorm(
    n = N,
    Means,
    CovMatrix,
    tol = 1e-6,
    empirical = FALSE,
    EISPACK = FALSE
  ))
  x <- rep("X",ncol(CovMatrix))
  names <- paste(x,1:ncol(CovMatrix),sep = "")
  names(simu_data) <- names
  return(simu_data)
}


