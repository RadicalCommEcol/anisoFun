
#' Random estimations of capital Omega
#'
#' Given an interaction matrix, this function generates a vector with 
#' estimations of capital Omega (i.e., the proportion of the feasible parameter 
#' space inside the unit sphere). Those estimations are computed via a 
#' quasi-Monte Carlo method for even relatively large communities.
#'
#' @param A_int a SxS interaction matrix, where S is the number of species.
#' @param replicates a number that specifies how many estimations of Omega
#'    will be calculated. By default, this parameter is set to 1,000 replicates.
#' @param use_chol_decomp boolean that specifies if the QR decomposition should
#' be used to compute Omega.
#' 
#' @import doParallel
#' @import foreach
#' @importFrom mvtnorm pmvnorm
#' @importFrom boot boot
#'
#' @return A vector whose length is equal to the parameter replicate, where each
#'  element represents an estimation of capital Omega.
#'
#' @examples
#' numCores <- detectCores()
#' numCores
#' registerDoParallel(2)
#' A_int <- -1*diag(c(3,5,7))
#' Omega_bootstrap(A_int, replicates = 1e3)
#' stopImplicitCluster()
#' 
#' @references \url{https://doi.org/10.1016/j.jtbi.2018.04.030}
#'
#' @export


Omega_bootstrap <- function(A_int, replicates = 1e3, use_chol_decomp = F) {
  
  # check the input parameters
  
  is.wholenumber <- function(x){
    integer_part <- round(x,0)
    decimal_part <- x - integer_part
    return(decimal_part == 0)
  }
  
  test_result_A_int <- check_A_int(A_int)
  test_replicates <- is.wholenumber(replicates)
  test_chol_decomp <- is.logical(use_chol_decomp)
  
  if((test_result_A_int[[1]]==T) & test_replicates & test_chol_decomp){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    if(use_chol_decomp == F){
      
      Sigma <- solve(t(A_int) %*% A_int,tol = 1e-17)
      
    }else{
      
      chol_A <- chol(t(A_int) %*% A_int)
      Sigma <- chol2inv(chol_A) 
      
    }
    
    S <- ncol(A_int)
    
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    
    res_list <- foreach (i = 1:replicates, .combine=c) %dopar% {
      
      d <- try(mvtnorm::pmvnorm(lower = rep(0, S),
                                upper = rep(Inf, S),
                                mean = rep(0, S), sigma = Sigma),silent = TRUE)
      out <- ifelse(class(d) == "try-error",0,d[1])
      
      
    }
    
    return(res_list)
    
    
  }else{
    
    if(test_result_A_int[[1]] != T){
      cat(test_result_A_int[[2]],"\n")
    }
    if(test_replicates != T){
      cat("replicates should be a whole number.","\n")
    }
    if(test_chol_decomp != T){
      cat("use_chol_decomp should be logical, that is, TRUE or FALSE.","\n")
    }
    
  }
  
}


