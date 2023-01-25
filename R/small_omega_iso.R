
#' Small Omega for the maximum isotropic cap
#'
#' Given an SxS interaction matrix, this function computes the value of small 
#' Omega for maximum isotropic cap (i.e., the S root of the maximum isotropic 
#' cap, Xi).
#'
#' @param A_int a SxS interaction matrix, where S is the number of species.
#' 
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#'
#' @return A number that shows the value of value of small Omega, the S root of 
#' the volume of the maximum isotropic cap, Xi.
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' small_omega_iso(A_int)
#'
#' @export

small_omega_iso <- function(A_int){
  
  test_result_A_int <- check_A_int(A_int)
  
  if(test_result_A_int[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    # coordinates of FD's Vertices
    vertices_unit_ball(A_int)
    
    # Get main curves data: vertices and normal vectors for their hyperplanes
    dimensions <- ncol(A_int)
    
    # Estimation of the incenter position and of the inradius for all the three 
    # vertices (of the simplicial cone) combinations
    
    incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
    
    isoprobability <- incenter_inradius_isoprob[[3]]
    
    return(isoprobability^(1/nrow(A_int)))
    
  }else{
    
    cat(test_result_A_int[[2]],"\n")
    
  }
  
}

