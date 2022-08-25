#source(system.file("analysis/stand_age/midnorth.stand_age.R", package="sibships"))

library(sibships)
library(raster)

load(system.file("analysis/data/midnorth.RData", package="sibships"))

landscape_covariates <- raster::stack(list("stand_age"=scale(stand_age)))

resistance_model <- function(raster_stack, parameters)
{
  # converts stack of rasters to a single output raster
  # representing a resistance surface
  stopifnot("stand_age" %in% names(raster_stack))
  stopifnot("theta" %in% names(parameters))
  stopifnot(length(parameters) == 1)
  resistance <- exp(raster_stack[["stand_age"]] * parameters["theta"])
  return(resistance)
}

parameter_grid <- as.matrix(expand.grid(
  "theta"=unique(c(
     seq(-1.5, -0.5, length.out=11), 
     seq(-0.5,  0.5, length.out=51),
     seq( 0.5,  1.5, length.out=11)
))))

#does model work? evaluate at first point in parameter grid, check for NAs, etc
plot(resistance_model(landscape_covariates, parameter_grid[1,]))

fit <- sibship_foraging_model(
  colony_count_at_traps, 
  floral_cover_at_traps, 
  trap_coordinates,
  landscape_covariates,
  resistance_model,
  parameter_grid,
  verbose=TRUE
)
save(fit, file="midnorth.stand_age.fitted.RData")

#simulate/refit at maximum likelihood estimates of the parameters
boot_at_mle <- parametric_bootstrap(fit, fit$mle, num_boot=100, verbose=TRUE, random_seed=1)

##simulate/refit at null model
null <- c("theta" = 0)
boot_at_null <- parametric_bootstrap(fit, null, num_boot=100, verbose=TRUE, random_seed=1)

save(boot_at_null, boot_at_mle, file="midnorth.stand_age.simulations.RData")

#figures
#plot_1d_likelihood_surface(fit, simulations=boot_at_mle) + 
#  ggtitle("Simulations from fitted model")
#plot_1d_likelihood_surface(fit, simulations=boot_at_null) + 
#  ggtitle("Simulations from null model")
#plot_1d_sampling_distributions(fit, parametric_bootstraps=boot_at_mle, null_simulations=boot_at_null) + 
#  ggtitle("Sampling distributions of the MLE")
