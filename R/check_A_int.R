
#' Internal, check the initial values of the interaction matrix
#'
#' @param A_int SxS interaction matrix, where S is the number of species.
#'
#' @return A list with 2 elements:
#' \itemize{
#'    \item The first element is a boolean that is equal to TRUE, if the
#'    matrix passes all tests.
#'    \item The second element is a string that contains an error or warning
#'    message.
#'  }
#' @noRd

check_A_int <- function(A_int){
  
  A_int_ok <- TRUE
  A_int_message <- ""
  
  # A_int should be a square matrix
  if(!all(class(A_int) %in% c("matrix","array"))){
    A_int_ok <- FALSE
    A_int_message <- "The interaction matrix should be matrix."
  }else if(!is.numeric(A_int)){
    A_int_ok <- FALSE
    A_int_message <- "The interaction matrix should be numeric."
  }else if(nrow(A_int)!=ncol(A_int)){
    A_int_ok <- FALSE
    A_int_message <- "The interaction matrix should be square."
  }else if(is.null(colnames(A_int)) & is.null(rownames(A_int))){
    A_int_message <- "Adding species labels to interaction matrix."
  }
  
  return(list(A_int_ok,A_int_message))
  
}