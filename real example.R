library(matlib) # to multiply matrices
library(tidyverse)
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(MultitrophicFun)
library(parallel)
library(foreach)
library(doParallel)
library(boot)


# Functions to run calculations about the isotropic area
pathnames <- list.files(pattern="[.]R$", path="R", full.names=TRUE);
sapply(pathnames, FUN=source)

# We create the metaweb for each year

matrix_entries_raw <- read_csv2("Data/caracoles_raw_data/alpha_heterogeneous_time.csv") %>%
  group_by(year,focal,neighbour) %>% summarise(magnitude = mean(magnitude, na.rm =T))


number_Omega_replicates <- 200
number_boot_replicates <- number_Omega_replicates

years_included <- matrix_entries_raw$year %>% unique()

prob_excl_caracoles <- NULL 

# Config parallelization
numCores <- detectCores()
registerDoParallel(numCores-1)


# Estimate prob_excl_caracoles
set.seed(123)

for(year_i in years_included){
  
  cat(year_i, "\n")
  
  A_int <- matrix_year_i(matrix_entries_raw, year_i)
  prob_excl_A_int <- prob_extinction_4_int_matrix(A_int, 
                                                  number_Omega_replicates, 
                                                  number_boot_replicates)
  
  prob_excl_A_int$year <- year_i
  prob_excl_caracoles <- bind_rows(prob_excl_caracoles, prob_excl_A_int)
  
}

stopImplicitCluster() #

number_plant_sp_year <- prob_excl_caracoles %>% group_by(year) %>% count() %>%
  rename(number_plant_sp = n)

# Add number of species per year, isoprobability of exclusion and mean exclusion ratio

caracoles_plot_year <- prob_excl_caracoles %>% 
  left_join(number_plant_sp_year, by = c("year")) %>%
  mutate(iso_prob_exclusion = 1/number_plant_sp,
         exclusion_ratio_mean = prob_excl_mean/iso_prob_exclusion,
         exclusion_ratio_upper = prob_excl_upperCI/iso_prob_exclusion,
         exclusion_ratio_lower = prob_excl_lowerCI/iso_prob_exclusion)

# Plot mean exclusion ratio

species_names <- c(
  `BEMA` = "Beta\nmacrocarpa",
  `CETE` = "Centaurium\ntenuiflorum",
  `HOMA` = "Hordeum\nmarinum",
  `LEMA` = "Leontodon\nmaroccanus",
  `PAIN` = "Parapholis\nincurva",
  `POMA` = "Polypogon\nmaritimus",
  `SASO` = "Salsola\nsoda",
  `SCLA` = "Scorzonera\nlaciniata"
)


ggplot(caracoles_plot_year %>% filter(species %in% c("BEMA","CETE", "HOMA",
                                                     "LEMA","PAIN","POMA",
                                                     "SASO")),
       aes(x = as.factor(year),y = exclusion_ratio_mean))+
  geom_errorbar(aes(ymin=exclusion_ratio_lower, ymax=exclusion_ratio_upper), width=.2)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.3)+
  stat_summary(fun = mean, geom="point", shape=23, size=3,fill="red")+
  geom_hline(yintercept = 1, color = "deepskyblue", linetype = "dashed", size = 1.3)+
  facet_wrap(~species, ncol = 4, labeller = as_labeller(species_names))+ #, scales = "free")+
  # scale_y_continuous(labels = scientific)+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  ylab("Exclusion ratio")+
  labs(x=NULL)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text = element_text(face = "italic"))
