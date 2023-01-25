
matrix_year_i <- function(matrix_entries_raw, year_i){
  
  data_year_i <- matrix_entries_raw %>% filter(year == year_i, !is.na(magnitude))
  
  
  if(nrow(data_year_i)>0){
    
    name_focals <- data_year_i$focal %>% unique() %>% sort()
    name_neighbors <- data_year_i$neighbour %>% unique() %>% sort()
    
    # if(length(name_focals)<length(name_neighbors)){
    #   name_plants <- name_focals
    # }else{
    #   name_plants <- name_neighbors
    # }
    
    name_plants <- intersect(name_focals,name_neighbors)
    
    plant_species <- length(name_plants)
    
    matrix_year_i_out <- matrix(rep(0,plant_species^2),
                                       nrow = plant_species, ncol = plant_species)
    
    row.names(matrix_year_i_out) <- name_plants
    colnames(matrix_year_i_out) <- name_plants
    
    for (focal_i in 1:plant_species) {
      
      for (neigh_j in 1:plant_species) {
        
        matrix_year_i_out[focal_i, neigh_j] <- data_year_i %>% 
          filter(focal == name_plants[focal_i],
                 neighbour == name_plants[neigh_j]) %>% ungroup() %>%
          select(magnitude) %>% as.numeric()
        
      }
      
    }
    
    matrix_year_i_out[is.na(matrix_year_i_out)] <- 0.0
    
    return(matrix_year_i_out)
    
  }else{
    
    matrix_year_i_out <- "Combination_absent"
    
    return(matrix_year_i_out)
  }
  
  
}
