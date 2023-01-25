generate_random_feasible_growth_rates <- function(A_int, random_points=1000){
  
  num_species <- ncol(A_int)
  inv_A_int <- inv(A_int)
  
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
  incenter <- incenter_inradius_isoprob[[1]]
  max_radius <- max_arc_distance_vertices_incenter(A_int, incenter)
  
  results <- data.frame(matrix(ncol = 2*ncol(A_int), nrow = random_points))
  r_names <- c(paste0("r_", 1:num_species),paste0("N_", 1:num_species))
  colnames(results) <- r_names
  
  number_feasible_points <- 0
  
  while(number_feasible_points < random_points){
    
    aux_growth_rate <- incenter+rnorm(num_species, mean = 0, sd = max_radius)
    growth_rate <- aux_growth_rate/sqrt(sum(aux_growth_rate*aux_growth_rate))
    abundance <- -inv_A_int %*% growth_rate %>% as.numeric()
    
    if(all(abundance>0)){
      number_feasible_points <- 1+number_feasible_points
      results[number_feasible_points,] <- c(growth_rate,abundance)
      
    }
    
  }
  
  return(results)
}


max_arc_distance_vertices_incenter <- function(A_int, incenter){
  
  number_species <- ncol(A_int)
  vertices <- vertices_unit_ball(A_int)
  
  vertices_growth_rates <- vertices[,(number_species+1):(2*number_species)]
  arc_distances <- NULL
  
  for(j in 1:number_species){
    
    arc_distances <- c(arc_distances,sum(as.numeric(vertices_growth_rates[j,])*incenter))
    
  }
  
  return(max(abs(arc_distances)))
}
