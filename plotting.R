
library(matlib) # to multiply matrices
library(tidyverse)
library(gg3D)
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(RColorBrewer)


# This function get the coordinates of npoints in a circle with radius "radius_circle", centered at
# the point C (located by "Center_vector") and in the plane that contain the vectors "unit_vector1",
# and "unit_vector2", where unit_vector1 and unit_vector2 are normal vectors
# (i.e., unit_vector1 . unit_vector2 = 0)

circle3DFun <- function(center_vector,radius_circle, unit_vector1,unit_vector2,npoints = 100){
  
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- as.numeric(center_vector[1]) + radius_circle * cos(tt) * unit_vector1[1] + 
    radius_circle * sin(tt) * unit_vector2[1]
  yy <- as.numeric(center_vector[2]) + radius_circle * cos(tt) * unit_vector1[2] + 
    radius_circle * sin(tt) * unit_vector2[2]
  zz <- as.numeric(center_vector[3]) + radius_circle * cos(tt) * unit_vector1[3] + 
    radius_circle * sin(tt) * unit_vector2[3]
  return(data.frame(x = xx, y = yy, z = zz))
}

# To plot circles
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# To plot the FD's limits
limits_feasibility_domain <- function(A_int, species_limiting, number_points_domain=1000){
  
  num_species <- 1:(ncol(A_int))
  
  r_values <- data.frame(matrix(ncol = (2*ncol(A_int)), nrow = number_points_domain))
  r_names <- c(paste0("N_", num_species),paste0("r_", num_species)) 
  colnames(r_values) <- r_names
  
  for(j in 1:ncol(A_int)){
    
    r_values[,j] <- runif(number_points_domain) # uniform on [-1, 1]
    
  }
  
  r_values[,species_limiting] <- 0
  
  for(i in 1:nrow(r_values)){
    
    N_column <- matrix(as.numeric(r_values[i,1:ncol(A_int)]),
                       ncol = 1, nrow = ncol(A_int) )
    
    r_column <- (-1)*(A_int %*% as.matrix(N_column))
    
    sq_sum <- sum(r_column*r_column)
    
    r_values[i,c((ncol(A_int)+1):(2*ncol(A_int)))] <- r_column/sqrt(sq_sum)
    
  }
  
  return(r_values)
  
}

# Create an interaction matrix
S = 3 # number of species

use_intencer <- TRUE #FALSE es un random point.


alphas <- runif(S*S)

set.seed(1234)
A_int <- -1*matrix(alphas, ncol=S, nrow=S)

# Get main curves data: vertices and normal vectors for their hyperplanes
dimensions <- ncol(A_int)
curves_data <- cone_vertices_director_vertices(A_int)

# Estimate the vector of initial intrinsic growth rates

if(use_intencer){
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
  center_unit_ball <- incenter_inradius_isoprob[[1]]
}else{
  
  # We choose a random point inside de feasibility domain (N_i^*>0) and estimate its growth rates
  set.seed(123)
  feasible_N <- runif(3)
  feasible_N_mat <- matrix(feasible_N,nrow = dimensions,ncol = 1)
  center_mat <- (-1)*A_int %*% feasible_N_mat
  
  center <- as.numeric(center_mat) 
  center_unit_ball <- center/sqrt(sum(center*center))
  
}



FD_point_closest_to_exclusion <- closest_to_exclusion(center_unit_ball,A_int)
min_great_circle_dist_to_limits <- acos(sum(center_unit_ball*as.numeric(FD_point_closest_to_exclusion[1,1:S])))


# Plot results---------
# 1) Data to plot the limits of the feasibility domain
r_N1 <- limits_feasibility_domain(A_int, species_limiting=1, number_points_domain=1000)
r_N2 <- limits_feasibility_domain(A_int, species_limiting=2, number_points_domain=1000)
r_N3 <- limits_feasibility_domain(A_int, species_limiting=3, number_points_domain=1000)

# 2) Data to plot the vertices of the simplicial cone, the incenter and the tangent points
important_points_feasibility_domain_aux <- vertices_unit_ball(A_int)
important_points_feasibility_domain_aux$color_sol <- "Vertex"
random_center_coordinates <- tibble(r_1 = center_unit_ball[1],
                                    r_2 = center_unit_ball[2],
                                    r_3 = center_unit_ball[3],
                                    color_sol = "Random Center")

important_points_feasibility_domain <- bind_rows(important_points_feasibility_domain_aux,
                                                 random_center_coordinates)
# 3) Data to plot a 2d Circle centered at the origin, i.e., the 2-D projection of the unit ball
circ_points <- circleFun(c(0,0),2,npoints = 100)

# 4) Data to plot a the incircle that is tangent to the feasibility domain
cos_theta <- cos(min_great_circle_dist_to_limits)
center_vector <- as.numeric(center_unit_ball*cos_theta)
radius_circle <- sin(min_great_circle_dist_to_limits)

point_1 <- c(0,0,cos_theta/center_vector[3])

unit_vector1 <- center_vector-point_1
unit_vector1 <- unit_vector1 / sqrt(sum(unit_vector1*unit_vector1))

# Cross product
unit_vector2 <- pracma::cross(unit_vector1,center_vector)
unit_vector2 <- unit_vector2 / sqrt(sum(unit_vector2*unit_vector2))

circ_3D_points <- circle3DFun(center_vector,radius_circle,
                              unit_vector1,unit_vector2,npoints = 100) %>%
  rename(r_1=x,r_2=y,r_3=z) %>% mutate(color_sol = "Tangent circle",
                                       unit_ball=r_1^2+r_2^2+r_3^2)

# 2D Plots
ggplot()+
  geom_path(data=circ_points,aes(x=x,y=y))+
  geom_path(data=circ_3D_points,aes(x=r_1,y=r_2))+
  geom_point(data=r_N1,aes(x=r_1,y=r_2),color="red")+
  geom_point(data=r_N2,aes(x=r_1,y=r_2),color="black")+
  geom_point(data=r_N3,aes(x=r_1,y=r_2),color="blue")+
  geom_point(data=important_points_feasibility_domain,aes(x=r_1,y=r_2),color="purple",size=2)+
  labs(x="Intrinsic growth rate sp 1", y="Intrinsic growth rate sp 2")+
  theme_bw()+
  theme(aspect.ratio=1)

ggplot()+
  geom_path(data=circ_points,aes(x=x,y=y))+
  geom_path(data=circ_3D_points,aes(x=r_1,y=r_3))+
  geom_point(data=r_N1,aes(x=r_1,y=r_3),color="red")+
  geom_point(data=r_N2,aes(x=r_1,y=r_3),color="black")+
  geom_point(data=r_N3,aes(x=r_1,y=r_3),color="blue")+
  geom_point(data=important_points_feasibility_domain,aes(x=r_1,y=r_3),color="purple",size=2)+
  labs(x="Intrinsic growth rate sp 1", y="Intrinsic growth rate sp 3")+
  theme_bw()+
  theme(aspect.ratio=1)

ggplot()+
  geom_path(data=circ_points,aes(x=x,y=y))+
  geom_path(data=circ_3D_points,aes(x=r_2,y=r_3))+
  geom_point(data=r_N1,aes(x=r_2,y=r_3),color="red")+
  geom_point(data=r_N2,aes(x=r_2,y=r_3),color="black")+
  geom_point(data=r_N3,aes(x=r_2,y=r_3),color="blue")+
  geom_point(data=important_points_feasibility_domain,aes(x=r_2,y=r_3),color="purple",size=2)+
  labs(x="Intrinsic growth rate sp 2", y="Intrinsic growth rate sp 3")+
  theme_bw()+
  theme(aspect.ratio=1)


# 3D Plots

r_N1$color_sol <- "N1=0"
r_N2$color_sol <- "N2=0"
r_N3$color_sol <- "N3=0"

feas_limits <- bind_rows(important_points_feasibility_domain,r_N1,r_N2,r_N3,circ_3D_points)

ggplot(data = feas_limits, aes(x=r_1,y=r_2, z=r_3,color=color_sol))+
  theme_void() +
  scale_color_brewer(palette = "Set1")+
  axes_3D() +
  stat_3D()+
  theme(aspect.ratio=1)+
  labs(color="Description",title="Feasibility domain (3 species)")
