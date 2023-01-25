#' Bootstrap estimation of the probabilities of exclusion and capital Omega
#'
#' Given an interaction matrix, this function generates a tibble with S rows,
#' where S is the number of species in the matrix. Each row contains the
#' following information for a given: replicates of its probability of exclusion
#' and the corresponding value of capital Omega for each replicate (i.e., the 
#' proportion of the feasible parameter #' space inside the unit sphere). Those 
#' estimations are computed via a quasi-Monte Carlo method and by 
#' bootstrapping the quasi-Mote Carlo results.
#'
#' @param A_int a SxS interaction matrix, where S is the number of species.
#' @param number_Omega_replicates a number that specify how many estimations 
#'    of Omega will be calculated by the quasi-Monte Carlo method. By default, 
#'    this parameter is set to 1,000 replicates.
#' @param number_boot_replicates a number that specifies how many bootstrap
#'    estimations of the probability of exclusion and Omega will be calculated
#'    for each species. By default, this parameter is set to 1,000 replicates.
#' @param use_chol_decomp boolean that specifies if the QR decomposition should
#' be used to compute Omega.
#' 
#' @import foreach
#' @import doParallel
#' @importFrom mvtnorm pmvnorm
#' @importFrom boot boot
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#'
#' @return A tibble with S rows and 2 x number_boot_replicates columns. Each row 
#' contains the estimated probabilities of exclusion and values of Omega for 
#' each species, when considering the intersection of the feasibility domain 
#' and the unit ball. Specifically, each row contains:
#'    \itemize{
#'    \item Resamples of the estimated probabilities of exclusion for species i 
#'    (prob_excl_rep_x), where x takes values from 1 to number_boot_replicates.
#'    \item The corresponding value of capital Omega for species i and 
#'    a given prob_excl_rep_x (Omega_rep_x), where x takes values from 1 to 
#'    number_boot_replicates.
#'    } 
#'    
#' @examples
#' numCores <- detectCores()
#' numCores
#' registerDoParallel(2)
#' A_int <- -1*diag(c(3,5,7))
#' boot_prob_excl_Omega_raw(A_int)
#' stopImplicitCluster()
#' 
#' @references \url{https://doi.org/10.1016/j.jtbi.2018.04.030}
#'
#' @export


boot_prob_excl_Omega_raw <- function(A_int, number_Omega_replicates = 1e3,
                                 number_boot_replicates = 1e3,
                                 use_chol_decomp = F){
  
  
  # check the input parameters
  
  is.wholenumber <- function(x){
    integer_part <- round(x,0)
    decimal_part <- x - integer_part
    return(decimal_part == 0)
  }
  
  test_result_A_int <- check_A_int(A_int)
  test_omega_rep <- is.wholenumber(number_Omega_replicates)
  test_boot_rep <- is.wholenumber(number_boot_replicates)
  test_chol_decomp <- is.logical(use_chol_decomp)
  
  if((test_result_A_int[[1]]==T) & test_omega_rep & test_boot_rep &
     test_chol_decomp){
    
    simple_mean <- function(x, indices){
      return(sum(x[indices])/length(indices))
    }
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    exclusion_probabilities_par_bootstrap <- function(A_int, I, replicates = 1e3,
                                                      use_chol_decomp = F){
      
      dimensions <- ncol(A_int)
      
      all_vertices_data_aux <- vertices_unit_ball(A_int)
      all_vertices_data <- all_vertices_data_aux[,c((dimensions+1):(2*dimensions))]
      
      feasibility_parts <- matrix(rep(0,nrow(A_int)*replicates),
                                  nrow = nrow(A_int), ncol = replicates)
      
      names_sp <- colnames(A_int)
      
      for (i in 1:nrow(A_int)) {
        
        cat(names_sp[i],"\n")
        
        A_int_mod <- A_int
        
        A_int_mod[,i] <- -I
        
        feasibility_parts[i, ] <- Omega_bootstrap(A_int_mod, replicates,
                                                  use_chol_decomp = F)
        
      }
      
      return(feasibility_parts)
      
    }
    
    
    dimensions <- ncol(A_int)
    
    incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
    
    I <- incenter_inradius_isoprob[[1]]
    
    name_colums <- c("species", paste0("prob_excl_rep_",1:number_boot_replicates),
                     paste0("Omega_rep_",1:number_boot_replicates))
    
    A_int_prob_exclusion_plants_mat <- matrix(rep(0, dimensions*(1 + 2*number_boot_replicates)),
                                              nrow = dimensions, 
                                              ncol = (1 + 2*number_boot_replicates))
    
    bootstrap_raw_data <-
      exclusion_probabilities_par_bootstrap(A_int, I,
                                            replicates = number_Omega_replicates,
                                            use_chol_decomp = F)
    
    for (sp_i in 1:dimensions) {
      
      bootmean_area_excl_sp_i <- boot::boot(data = bootstrap_raw_data[sp_i,], 
                                            statistic = simple_mean, 
                                            R = number_boot_replicates)
      
      A_int_prob_exclusion_plants_mat[sp_i, 2:(1+number_boot_replicates)] <- 
        bootmean_area_excl_sp_i[["t"]] 
      
    }
    
    Omega_values <- colSums(A_int_prob_exclusion_plants_mat[, c(2:(1+number_boot_replicates))])
    
    A_int_prob_exclusion_plants_mat[, c((2+number_boot_replicates):(1 + 2*number_boot_replicates))] <-
      matrix(rep(Omega_values,dimensions), nrow = number_boot_replicates) %>% t()
    
    for (sp_i in 1:nrow(A_int_prob_exclusion_plants_mat)) {
      A_int_prob_exclusion_plants_mat[sp_i, c(2:(1+number_boot_replicates))] <-
        A_int_prob_exclusion_plants_mat[sp_i, c(2:(1+number_boot_replicates))]/
        A_int_prob_exclusion_plants_mat[sp_i, c((2+number_boot_replicates):(1 + 2*number_boot_replicates))]
    }
    
    A_int_prob_exclusion_plants_aux <- tibble::as_tibble(data.frame(A_int_prob_exclusion_plants_mat))
    colnames(A_int_prob_exclusion_plants_aux) <- name_colums
    
    A_int_prob_exclusion_plants_aux$species = colnames(A_int)
    
    return(A_int_prob_exclusion_plants_aux)
    
  }else{
    
    if(test_result_A_int[[1]] != T){
      cat(test_result_A_int[[2]],"\n")
    }
    if(test_omega_rep != T){
      cat("number_Omega_replicates should be a whole number.","\n")
    }
    if(test_boot_rep != T){
      cat("number_boot_replicates should be a whole number.","\n")
    }
    if(test_chol_decomp != T){
      cat("use_chol_decomp should be logical, that is, TRUE or FALSE.","\n")
    }
  }
  
}
