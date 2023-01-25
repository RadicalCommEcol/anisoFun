context("parameter fitting functions")

# set data ----------------------------------------------------------------

S = 3
A_int_symm <- (-1)*diag(rep(1,S))
incenter_inradius_isoprob_sym <- 
  incenter_inradius_isoprob_calculation(A_int_symm)
I_sym <- incenter_inradius_isoprob_sym[[1]]
Xi_sym <- incenter_inradius_isoprob_sym[[3]]
  
center_rand_aux <- runif(S)
center_rand <- center_rand_aux/sqrt(sum(center_rand_aux^2))
  
# Interaction matrices with the same value of capital Omega


A_int_1 <- matrix(rep(0,9), ncol=3)
A_int_1[,1] <- -c(1, 0, 0)

A_int_2 <- A_int_1

A_int_1[,2] <- -c(0.3406884710289558, 0.9328944522580684, 0.11678744219579273)
A_int_1[,3] <- -c(0.3406884710289558, 0.12410827817548943, 0.9319487652208257)

A_int_2[,2] <- -c(0.3995691276168492, 0.7476891904854146, -0.5303823023903146)
A_int_2[,3] <- -c(0.39956912762129915, -0.19648768029759267, 0.8953977349474517)

incenter_inradius_isoprob_mat_1 <- 
  incenter_inradius_isoprob_calculation(A_int_1)
I_1 <- incenter_inradius_isoprob_mat_1[[1]]
theta_1 <- incenter_inradius_isoprob_mat_1[[2]]
Xi_1 <- incenter_inradius_isoprob_mat_1[[3]]

incenter_inradius_isoprob_mat_2 <- 
  incenter_inradius_isoprob_calculation(A_int_2)
I_2 <- incenter_inradius_isoprob_mat_2[[1]]
theta_2 <- incenter_inradius_isoprob_mat_2[[2]]
Xi_2 <- incenter_inradius_isoprob_mat_2[[3]]


# when run with nonsensical parameters, the function returns NULL

test_that("invalid arguments return NULL",{
  
  wrong_A_int_1 <- c(2,"b",1,"x")
  
  wrong_anisotropy_index_1 <- anisotropy_index(wrong_A_int_1)
  wrong_mean_pairwise_sim_1 <- mean_pairwise_similarity(wrong_A_int_1)
  
  wrong_anisotropy_metrics_1 <- anisotropy_metrics_4_int_matrix(wrong_A_int_1)
  wrong_boot_prob_excl_Omega_1 <- boot_prob_excl_Omega_raw(wrong_A_int_1)
  wrong_cone_vertices_director_vertices_1 <- cone_vertices_director_vertices(wrong_A_int_1)
  wrong_incenter_inradius_isoprob_1 <- incenter_inradius_isoprob_calculation(wrong_A_int_1)
  wrong_Omega_bootstrap_1 <- Omega_bootstrap(wrong_A_int_1)
  wrong_prob_extinction_1 <- prob_extinction_4_int_matrix(wrong_A_int_1)
  wrong_small_omega_iso_1 <- small_omega_iso(wrong_A_int_1)
  wrong_vertices_1 <- vertices_unit_ball(wrong_A_int_1)
  wrong_closest_to_excl_1 <- closest_to_exclusion(center_rand, wrong_A_int_1)
  wrong_furthest_from_excl_1 <- furthest_from_exclusion(center_rand, wrong_A_int_1)
  

  wrong_A_int_2 <- matrix(c(1,2,3,"x"),ncol = 2)
  
  wrong_anisotropy_index_2 <- anisotropy_index(wrong_A_int_2)
  wrong_mean_pairwise_sim_2 <- mean_pairwise_similarity(wrong_A_int_2)
  wrong_anisotropy_metrics_2 <- anisotropy_metrics_4_int_matrix(wrong_A_int_2)
  wrong_boot_prob_excl_Omega_2 <- boot_prob_excl_Omega_raw(wrong_A_int_2)
  wrong_cone_vertices_director_vertices_2 <- cone_vertices_director_vertices(wrong_A_int_2)
  wrong_incenter_inradius_isoprob_2 <- incenter_inradius_isoprob_calculation(wrong_A_int_2)
  wrong_closest_to_excl_2 <- closest_to_exclusion(center_rand, wrong_A_int_2)
  wrong_furthest_from_excl_2 <- furthest_from_exclusion(center_rand, wrong_A_int_2)
  wrong_Omega_bootstrap_2 <- Omega_bootstrap(wrong_A_int_2)
  wrong_prob_extinction_2 <- prob_extinction_4_int_matrix(wrong_A_int_2)
  wrong_small_omega_iso_2 <- small_omega_iso(wrong_A_int_2)
  wrong_vertices_2 <- vertices_unit_ball(wrong_A_int_2)
  
  
  wrong_probability_1 <- c(-.25,.75)
  wrong_anisotropy_index_3 <- anisotropy_index(wrong_probability_1)
  wrong_mean_pairwise_sim_3 <- mean_pairwise_similarity(wrong_probability_1)
  
  wrong_probability_2 <- c(.25,1)
  wrong_anisotropy_index_4 <- anisotropy_index(wrong_probability_2)
  wrong_mean_pairwise_sim_4 <- mean_pairwise_similarity(wrong_probability_2)
  

  expect_null(wrong_anisotropy_index_1)
  expect_null(wrong_mean_pairwise_sim_1)
  expect_null(wrong_closest_to_excl_1)
  expect_null(wrong_furthest_from_excl_1)
  expect_null(wrong_anisotropy_metrics_1)
  expect_null(wrong_boot_prob_excl_Omega_1)
  expect_null(wrong_cone_vertices_director_vertices_1)
  expect_null(wrong_incenter_inradius_isoprob_1)
  expect_null(wrong_Omega_bootstrap_1)
  expect_null(wrong_prob_extinction_1)
  expect_null(wrong_small_omega_iso_1)
  expect_null(wrong_vertices_1)
  
  expect_null(wrong_anisotropy_index_2)
  expect_null(wrong_mean_pairwise_sim_2)
  expect_null(wrong_closest_to_excl_2)
  expect_null(wrong_furthest_from_excl_2)
  expect_null(wrong_anisotropy_metrics_2)
  expect_null(wrong_boot_prob_excl_Omega_2)
  expect_null(wrong_cone_vertices_director_vertices_2)
  expect_null(wrong_incenter_inradius_isoprob_2)
  expect_null(wrong_Omega_bootstrap_2)
  expect_null(wrong_prob_extinction_2)
  expect_null(wrong_small_omega_iso_2)
  expect_null(wrong_vertices_2)
  
  
  expect_null(wrong_anisotropy_index_3)
  expect_null(wrong_mean_pairwise_sim_3)
  expect_null(wrong_anisotropy_index_4)
  expect_null(wrong_mean_pairwise_sim_4)

  
})

test_that("isotropic communities return expected values",{
  
  iso_prob <- rep(1/S,S)
  
  iso_anisotropy_index <- anisotropy_index(iso_prob)
  iso_mean_pairwise_sim <- mean_pairwise_similarity(iso_prob)
  iso_small_omega_iso <- small_omega_iso(A_int_symm)
  iso_closest_to_exclusion <- closest_to_exclusion(I_sym, A_int_symm)
  iso_furthest_from_exclusion <- furthest_from_exclusion(I_sym, A_int_symm)
  
  
  expect_equal(iso_anisotropy_index,1)
  expect_equal(iso_mean_pairwise_sim,1)
  expect_equal(iso_small_omega_iso,Xi_sym^(1/S))
  expect_equal(nrow(iso_closest_to_exclusion),S)
  expect_equal(length(iso_furthest_from_exclusion),S)
  
  
  skip_on_cran()
  
  registerDoParallel(2)
  iso_prob_extinction <- prob_extinction_4_int_matrix(A_int_symm)
  iso_Omega_bootstrap <- Omega_bootstrap(A_int_symm)
  stopImplicitCluster()
  
  expect_equal(round(dplyr::pull(iso_prob_extinction[,2]),3),round(iso_prob,3))
  expect_equal(round(dplyr::pull(iso_prob_extinction[,3]),3),round(iso_prob,3))
  expect_equal(round(dplyr::pull(iso_prob_extinction[,4]),3),round(iso_prob,3))
  expect_equal(length(iso_Omega_bootstrap),1e3)
  
  
})

test_that("given two comunities with the same Omega, the more isotropic has larger theta and Xi",{
  
  
  expect_lt(theta_2, theta_1)
  expect_lt(Xi_2, Xi_1)
  
  small_omega_iso_1 <- small_omega_iso(A_int_1)
  small_omega_iso_2 <- small_omega_iso(A_int_2)
  
  expect_lt(small_omega_iso_2, small_omega_iso_1)
  
  skip_on_cran()
  registerDoParallel(2)
  anisotropy_metrics_1 <- anisotropy_metrics_4_int_matrix(A_int_1)
  anisotropy_metrics_2 <- anisotropy_metrics_4_int_matrix(A_int_2)
  stopImplicitCluster()
  
  expect_equal(round(dplyr::pull(anisotropy_metrics_2[1,1:3]),3),round(dplyr::pull(anisotropy_metrics_1[1,1:3]),3)) # Same capital Omega
  expect_lt(round(dplyr::pull(anisotropy_metrics_2[1,4:6]),3),round(dplyr::pull(anisotropy_metrics_1[1,4:6]),3)) # Different shape
  
})