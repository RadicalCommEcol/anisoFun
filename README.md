[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8086505.svg)](https://doi.org/10.5281/zenodo.8086505)

# anisoFun 1.0.0

anisoFun provides a complete toolbox for calculating i) the maximum isotropic cap, ii) average exclusion probabilities, iii) anisotropic index, and iv) values of (small and capital) Omega for a given multispecies community, whose population dynamics can be approximated by a Lotka-Volterra (LV) model. The functions report 95% confidence intervals for those features that are estimated via quasi-Monte Carlo methods.

### Installation

The development version can be installed via the `remotes` package:

```R
install.packages("devtools")
require(devtools)
install_github("RadicalCommEcol/anisoFun")
library(anisoFun)
```

### Key features

The package has several key functions:

- `anisotropy_metrics_4_int_matrix()`
- `prob_extinction_4_int_matrix()`
- `incenter_inradius_isoprob_calculation()` 
- `closest_to_exclusion()`

Firstly, we create an interaction matrix for a 3 species community

```R
A_int <- -1*diag(c(3,5,7))
colnames(A_int) <- paste0("sp_",1:nrow(A_int))
rownames(A_int) <- paste0("sp_",1:nrow(A_int))
```
We estimate the following parameters:
* `I`: the incenter of the feasibility domain.
* `theta`: the great circle distance between the incerter and the feasibility domain's border.
* `Xi`: the maximum isotropic cap that can be placed in the feasibility domain.

```R
incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
I <- incenter_inradius_isoprob[[1]]
theta <- incenter_inradius_isoprob[[2]]
Xi <- incenter_inradius_isoprob[[3]]
```

In case you want to calculate the value of small `omega_iso` for maximum isotropic cap (i.e., the S root of the maximum isotropic cap, Xi), just run:

```R
omega_iso <- small_omega_iso(A_int)
```

For the above community, we can also estimate:
* The mean value of (small) `omega`, the `S` root of capital `Omega` (the proportion of the feasible parameter space inside the unit sphere), where `S` is the number of species in the community. We also provide the lower and upper bounds of the 95% confidence interval for small `omega`.
* The mean value of the anisotropy index `J`, along with the lower and upper bounds of the corresponding 95% confidence interval.

To reduce calculation time, we rely on parallel computation. In next example, we show how to check the number of cores in our computer (in our case, `numCores` = 12 cores) and how to select -for instance- `numCores` minus 4 cores.
```R
# Select the number of cores to perform a parallel computation and reduce calculation time
numCores <- detectCores()
numCores_selected <- numCores-4
registerDoParallel(numCores_selected)
```
If you skip the number of cores selection, one core will be used by default.

Once you selected the number of cores, to estimate (small) `omega` and the anisotropy index `J` and their confidence intervals, just run:

```R
anisotropy_metrics_A_int <- anisotropy_metrics_4_int_matrix(A_int)
anisotropy_metrics_A_int
```
To obtain the species' average probabilities of exclusion, execute the following command:

```R
prob_extinction_sp <- prob_extinction_4_int_matrix(A_int)
prob_extinction_sp
```

Finally, given an initial vector of intrinsic growing rates, denoted by center, we can compute the intrinsic growing rates and abundances of the species that is/are closest to exclusion, by running:

```R
center <- I + runif(nrow(A_int))
closect_sp_to_exclusion <- closest_to_exclusion(center, A_int)
closect_sp_to_exclusion
```





