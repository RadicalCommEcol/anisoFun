# load necessary packages

library(mvtnorm)
library(mgcv)
if(!require(ggtern)) {install.packages("ggtern"); library(ggtern)}
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}



first_species_excluded_ODE <- function(A_int, growth_rates_vector,initial_abundances, delta_t = 0.1, tmax = 1000) {
  
  number_species <- ncol(A_int)
  species_names <- colnames(A_int)
  closest_to_exclusion_estimation <- closest_to_exclusion(growth_rates_vector,A_int)
  perturbed_growth_rates_vector <- closest_to_exclusion_estimation[1:number_species] %>%
    as.numeric()
  
  alpha <- A_int
  dimnames(alpha) <- NULL
  
  change_growth_vector <- perturbed_growth_rates_vector-growth_rates_vector
  
  time_step <- seq(0,tmax,by=delta_t)
  r <- growth_rates_vector+1*change_growth_vector
  N0 <- initial_abundances
  parms <- list(r=r, alpha = alpha) #ODE
  model <- function(t,N,parms){ dN <- N * (parms$r + parms$alpha %*% N); list(dN)}
  sol <- ode(N0,time_step,model,parms, method = 'bdf', rtol = 1e-15, atol = 1e-15, maxsteps = 10000)# ode(N0,time_step,model,parms)
  # plot(sol)
  species_excluded <- which(min(sol[nrow(sol),2:(1+number_species)])==sol[nrow(sol),2:(1+number_species)]) %>% unname()
  
  return(species_names[species_excluded])
}
