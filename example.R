#load several functions
library(matlib) # to multiply matrices
library(tidyverse)
library(dplyr)
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(parallel)
library(foreach)
library(doParallel)
library(boot)

#This is to load all functions from AnisoFUN package
pathnames <- list.files(pattern="[.]R$", path="R", full.names=TRUE);
sapply(pathnames, FUN=source)

#The first step is to provide a matrix of species interactions

# Define symmetric interaction matrix--------------------------------------------
A_int <- matrix(rep(0,9), ncol=3)
A_int[,1] <- -c(1, 0, 0)
A_int[,2] <- -c(0.3406884710289558, 0.9328944522580684, 0.11678744219579273)
A_int[,3] <- -c(0.3406884710289558, 0.12410827817548943, 0.9319487652208257)

species_names <- paste0("sp_",1:ncol(A_int))
colnames(A_int) <- species_names
rownames(A_int) <- species_names
A_int

numCores <- detectCores()
cl <- makeCluster(numCores -1)
registerDoParallel(cl) # register the cluster for using foreach

prob_excl_AnisoFun_symmetric <- prob_extinction_4_int_matrix(A_int)
prob_excl_AnisoFun_symmetric

incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
incenter_inradius_isoprob

I <- incenter_inradius_isoprob[[1]] #the incenter of the feasibility domain.
theta <- incenter_inradius_isoprob[[2]] #the great circle distance between the incerter and the feasibility domain's border.
Xi <- incenter_inradius_isoprob[[3]] #the maximum isotropic cap that can be placed in the feasibility domain.


omega_iso <- small_omega_iso(A_int) # the S root of the maximum isotropic cap, Xi
omega_iso

numCores <- detectCores()
numCores_selected <- numCores-4
registerDoParallel(numCores_selected)


anisotropy_metrics_A_int <- anisotropy_metrics_4_int_matrix(A_int)
anisotropy_metrics_A_int

#Here
# small omega is the proportion of the feasible parameter space inside the unit sphere
# The mean value of the anisotropy index J,

prob_extinction_sp <- prob_extinction_4_int_matrix(A_int)
prob_extinction_sp

source("R/anisotropy_index.R")
probability_vector <- as.numeric(as.matrix(prob_extinction_sp[2]))
anisotropy_index<-anisotropy_index(probability_vector)
anisotropy_index

source("R/mean_pairwise_similarity.R")
mean_pairwise_similarity<-mean_pairwise_similarity(probability_vector)
mean_pairwise_similarity

center <- I + runif(nrow(A_int))
closest_sp_to_exclusion <- closest_to_exclusion(center, A_int)
closest_sp_to_exclusion

# Exclusion prob from ODE to compare results and speed. Careful, it takes time!
number_random_points <- 5e4

number_species <- ncol(A_int)

set.seed(1234)
random_feasible_growth_rates <- generate_random_feasible_growth_rates(A_int, random_points=number_random_points)

first_sp_excluded_results_symmetric <- foreach (i=1:nrow(random_feasible_growth_rates), .combine=c,  
                                                .packages = c("tidyverse","matlib",
                                                              "zipfR","pracma","deSolve")) %dopar% {
                                                                growth_rates_vector <- random_feasible_growth_rates[i,1:number_species] %>% as.numeric() %>% unname()
                                                                initial_abundances <- random_feasible_growth_rates[i,(number_species+1):(2*number_species)] %>% as.numeric() %>% unname()
                                                                species_excluded_i <- first_species_excluded_ODE(A_int, growth_rates_vector,initial_abundances, delta_t = 0.1, tmax = 1000)
                                                              }


random_feasible_growth_rates_results <- random_feasible_growth_rates
random_feasible_growth_rates_results$LV_first_sp_excluded <- first_sp_excluded_results_symmetric

prob_excl_ODE_symmetric <- table(first_sp_excluded_results_symmetric)/number_random_points
prob_excl_ODE_symmetric

stopImplicitCluster()




# Define asymmetric interaction matrix--------------------------------------------
A_int2 <- matrix(rep(0,9), ncol=3)
A_int2[,1] <- -c(1,0,0)
A_int2[,2] <- -c(0.3995691276168492, 0.7476891904854146, -0.5303823023903146)
A_int2[,3] <- -c(0.39956912762129915, -0.19648768029759267, 0.8953977349474517)

species_names <- paste0("sp_",1:ncol(A_int2))
colnames(A_int2) <- species_names
rownames(A_int2) <- species_names

# Exclusion prob from Anisofun
prob_excl_AnisoFun_asymmetric <- prob_extinction_4_int_matrix(A_int2)
prob_excl_AnisoFun_asymmetric

incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int2)
incenter_inradius_isoprob

I <- incenter_inradius_isoprob[[1]] #the incenter of the feasibility domain.
theta <- incenter_inradius_isoprob[[2]] #the great circle distance between the incerter and the feasibility domain's border.
Xi <- incenter_inradius_isoprob[[3]] #the maximum isotropic cap that can be placed in the feasibility domain.


omega_iso <- small_omega_iso(A_int2) # the S root of the maximum isotropic cap, Xi
omega_iso

numCores <- detectCores()
numCores_selected <- numCores-4
registerDoParallel(numCores_selected)


anisotropy_metrics_A_int <- anisotropy_metrics_4_int_matrix(A_int2)
anisotropy_metrics_A_int

#Here
# small omega is the proportion of the feasible parameter space inside the unit sphere
# The mean value of the anisotropy index J,

prob_extinction_sp <- prob_extinction_4_int_matrix(A_int2)
prob_extinction_sp

source("R/anisotropy_index.R")
probability_vector <- as.numeric(as.matrix(prob_extinction_sp[2]))
anisotropy_index<-anisotropy_index(probability_vector)
anisotropy_index

mean_pairwise_similarity<-mean_pairwise_similarity(probability_vector)
mean_pairwise_similarity

center <- I + runif(nrow(A_int2))
closest_sp_to_exclusion <- closest_to_exclusion(center, A_int2)
closest_sp_to_exclusion


# Increase the number of species to 6 and 12 species
#6 species
A_int2 <- matrix(rnorm(36),nrow=6)
diag(A_int2)<- -4
colnames(A_int2) <- paste0("sp_",1:nrow(A_int2))
rownames(A_int2) <- paste0("sp_",1:nrow(A_int2))

incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int2)
incenter_inradius_isoprob

I <- incenter_inradius_isoprob[[1]] #the incenter of the feasibility domain.
theta <- incenter_inradius_isoprob[[2]] #the great circle distance between the incerter and the feasibility domain's border.
Xi <- incenter_inradius_isoprob[[3]] #the maximum isotropic cap that can be placed in the feasibility domain.


omega_iso <- small_omega_iso(A_int2) # the S root of the maximum isotropic cap, Xi
omega_iso

numCores <- detectCores()
numCores_selected <- numCores-4
registerDoParallel(numCores_selected)


anisotropy_metrics_A_int <- anisotropy_metrics_4_int_matrix(A_int2)
anisotropy_metrics_A_int

#Here
# small omega is the proportion of the feasible parameter space inside the unit sphere
# The mean value of the anisotropy index J,

prob_extinction_sp <- prob_extinction_4_int_matrix(A_int2)
prob_extinction_sp

source("R/anisotropy_index.R")
probability_vector <- as.numeric(as.matrix(prob_extinction_sp[2]))
anisotropy_index<-anisotropy_index(probability_vector)
anisotropy_index

mean_pairwise_similarity<-mean_pairwise_similarity(probability_vector)
mean_pairwise_similarity

center <- I + runif(nrow(A_int2))
closest_sp_to_exclusion <- closest_to_exclusion(center, A_int2)
closest_sp_to_exclusion

#12 species
A_int3 <- matrix(rnorm(144),nrow=12)
diag(A_int3)<- -3
colnames(A_int3) <- paste0("sp_",1:nrow(A_int3))
rownames(A_int3) <- paste0("sp_",1:nrow(A_int3))

incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int3)
incenter_inradius_isoprob
I <- incenter_inradius_isoprob[[1]]
theta <- incenter_inradius_isoprob[[2]]
Xi <- incenter_inradius_isoprob[[3]]


omega_iso <- small_omega_iso(A_int3)
omega_iso

numCores <- detectCores()
numCores_selected <- numCores-4
registerDoParallel(numCores_selected)

source("R/anisotropy_index.R")
anisotropy_metrics_A_int <- anisotropy_metrics_4_int_matrix(A_int3)
anisotropy_metrics_A_int

prob_extinction_sp <- prob_extinction_4_int_matrix(A_int3)
prob_extinction_sp

probability_vector <- as.numeric(as.matrix(prob_extinction_sp[2]))
anisotropy_index<-anisotropy_index(probability_vector)
anisotropy_index

source("R/mean_pairwise_similarity.R")
mean_pairwise_similarity<-mean_pairwise_similarity(probability_vector)
mean_pairwise_similarity

center <- I + runif(nrow(A_int3))
closest_sp_to_exclusion <- closest_to_exclusion(center, A_int3)
closest_sp_to_exclusion
###
