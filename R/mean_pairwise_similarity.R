
#' Mean pairwise similarity
#'
#' Given a vector whose elements are the average probabilities of extinction of
#' the species in a given community, this function calculates the mean pairwise
#' similarity of such probabilities.
#' 
#' @param probability_vector a vector with length S whose elements are the 
#' average probabilities of extinction of the species in a given community, 
#' where S is the number of species.
#'
#' @return A non-negative number between 0 and 1 (isotropic community).
#'    
#' @examples
#' probability_vector <- c(.25,.15,.60)
#' mean_pairwise_similarity(probability_vector)
#'
#' @export

mean_pairwise_similarity <- function(probability_vector){
  
  # check the input parameter: probability_vector
  
  if(!all(class(probability_vector) %in% c("numeric"))){
    cat("Probabilities should be numbers.\n")
  }else if(all(class(probability_vector) %in% c("matrix","array"))){
    cat("The input parameter should be vector.\n")
  }else if(!all(probability_vector >=0)){
    cat("Probabilities should be non-negaive numbers.\n")
  }else if(round(sum(probability_vector),2) != 1){
    cat("The sum of all probabilities should be equal to one.\n")
  }else{
    
    # If the vector with probabilities is OK, then the index is calculated.
    
    number_species <- length(probability_vector)
    
    complement <- 0
    
    for(i in 1:(number_species-1)){
      
      for(j in (i+1):number_species){
        
        complement <- complement + abs(probability_vector[i]-probability_vector[j])
        
      }
      
    }
    
    return(1-complement/number_species)
    
  }
  
}
