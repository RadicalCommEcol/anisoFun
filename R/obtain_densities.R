#' Density of the species at a given intrinsic growth rates vector inside or outside the feasibility domain
#'
#' Given a vector r_norm2 (so norm2(r_norm2) = 1), computes the densities of the species at equilibrium
#'
#' @param r_norm2 an S vector that contains an intrinsic growth rate vector with sqrt(sum(r_norm2*r_norm2)) = 1
#' @param A_int SxS interaction matrix, where S is the number of species.
#'
#' @import quadprog
#' @return A matrix with 1 row and S columns. Each column contains the density for each species
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' r_norm2 <- matrix(c(1,1,1)/sqrt(3),ncol = 1)
#' obtain_densities(r_norm2, A_int) #output should be the (1/3,1/3,1/3)
#'
#' @export


obtain_densities <- function(r_norm2, A_int){
  
  test_result_A_int <- check_A_int(A_int)
  test_result_r_norm2 <- check_r_norm2(r_norm2, A_int)
  
  if(test_result_A_int[[1]]==T & test_result_r_norm2[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    if(test_result_r_norm2[[2]] == "Turning the initial growth rates into a unit vector."){
      
      r_norm2 <- r_norm2/sqrt(sum(r_norm2*r_norm2))
      
      cat(test_result_A_int[[2]],"\n")
      
    }

  # Using quadratic programming
  G <- t(rbind(-A_int,diag(length(r_norm2))))
  h <- matrix(c(r_norm2,rep(0,length(r_norm2))))
  sol <- solve.QP(-t(A_int+t(A_int)), r_norm2, G, bvec = h)
  N_star <- round(sol$solution, digits=12)
  
  #Project N_star to the simplex
  N_star <- N_star / sum(N_star)

  return(N_star)

  }else{
    
    if(test_result_A_int[[1]] != T){
      cat(test_result_A_int[[2]],"\n")
    }
    if(test_result_r_norm2[[1]] != T){
      cat(test_result_r_norm2[[2]],"\n")
    }
  }   
}
    
