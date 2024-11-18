#' Density of the species at a given intrinsic growth rates vector inside the feasibility domain
#'
#' Given
#'
#' @param center an S vector that contains an intrinsic growth rate vector 
#' @param A_int SxS interaction matrix, where S is the number of species.
#'
#' @return A matrix with 1 row and S columns. Each column contains the density for each species
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' obtain_densities(A_int,center)
#'
#' @export

obtain_densities <- function(center, A_int){
  
  test_result_A_int <- check_A_int(A_int)
  test_result_center <- check_center(center, A_int)
  
  if(test_result_A_int[[1]]==T & test_result_center[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    if(test_result_center[[2]] == "Turning the initial growth rates into a unit vector."){
      
      center <- center/sqrt(sum(center*center))
      
      cat(test_result_A_int[[2]],"\n")
      
    }

  } else {
    
    if(test_result_center[[2]] == "The vector that contains the initial growth rates should be in the feasibility domain.")){
      error_messages <- paste0(test_result_A_int[[2]],"\n",
                               test_result_center[[2]],"\n")
      
    }
    cat(error_messages)
  }
  
dimensions <- ncol(A_int)    
num_species <- 1:dimensions
N_star <- matrix(0,nrow = 1, ncol = dimensions)
sp_names <- paste0("N_", num_species)
colnames(N_star) <- sp_names
    
for(i in 1:dimensions){
  #Turning each column of A into a unit vector
  sup_A_col_i <- A_int[,i]
  sup_A_col_i <- sup_A_col_i / sqrt(sum(sup_A_col_i*sup_A_col_i))
  
  N_star[i] <- - center[i] / sup_A_col_i[i]
  

  }

  return(N_star)

}
    
