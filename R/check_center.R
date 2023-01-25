#' Internal, check the initial values of the interaction matrix
#'
#' @param center a vector with S elements, where S is the number of species.
#' @param A_int SxS interaction matrix, where S is the number of species.
#' 
#' @return A list with 2 elements:
#' \itemize{
#'    \item The first element is a boolean that is equal to TRUE, if the
#'    vector passes all tests.
#'    \item The second element is a string that contains an error or warning
#'    message.
#'  }
#' @noRd

check_center <- function(center, A_int){
  
  test_result_A_int <- check_A_int(A_int)
  
  if(test_result_A_int[[1]]==T){
    
    center_ok <- TRUE
    center_message <- ""
    
    populations_center <- (-1)*matlib::inv(A_int) %*% center
    module_center <- sqrt(sum(center*center))
    
    if(!is.numeric(center)){
      center_ok <- FALSE
      center_message <- "The vector that contains the initial growth rates should be numeric."
    }else if(nrow(A_int)!=length(center)){
      center_ok <- FALSE
      center_message <- "The vector that contains the initial growth rates and the interaction matrix should have the same number of species."
    }else if(!all(populations_center > 0)){
      center_ok <- FALSE
      center_message <- "The vector that contains the initial growth rates should be in the feasibility domain."
    }else if(module_center != 1){
      center_message <- "Turning the initial growth rates into a unit vector."}
    
    return(list(center_ok,center_message))
    
  }else{
    
    center_ok <- FALSE
    center_message <- "Invalid matrix does not allow to test the center"
    
    return(list(center_ok,center_message))
    
  }
  
  
 }


