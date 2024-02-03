#source(system.file("analysis/stand_age/north.stand_age.R", package="sibships"))

library(sibships)
library(raster)

load(system.file("analysis/data/midnorth.RData", package="sibships"))

roads <- stand_age*0 + sample(0:1, ncell(stand_age), replace=TRUE)

landscape_covariates <- raster::stack(
  list("roads"=roads, "age"=scale(stand_age))
)

resistance_model <- function(raster_stack, parameters)
{
  # converts stack of rasters to a single output raster
  # representing a resistance surface
  stopifnot("age" %in% names(raster_stack))
  stopifnot("roads" %in% names(raster_stack))
  stopifnot("theta_age" %in% names(parameters))
  stopifnot("theta_roads" %in% names(parameters))
  stopifnot(length(parameters) == 2)
  resistance <- exp(
    raster_stack[["age"]] * parameters["theta_age"] +
    raster_stack[["roads"]] * parameters["theta_roads"]
  )
  return(resistance)
}

parameter_grid_full <- as.matrix(expand.grid(
  "theta_age"=seq(-0.5,  0.5, length.out=3),
  "theta_roads"=seq(-0.5, 0.5, length.out=3)
))

parameter_grid_null <- as.matrix(expand.grid(
  "theta_age"=seq(-0.5,  0.5, length.out=11),
  "theta_roads"=0
))

#does model work? evaluate at first point in parameter grid, check for NAs, etc
plot(resistance_model(landscape_covariates, parameter_grid_null[1,]))
plot(resistance_model(landscape_covariates, parameter_grid_full[1,]))

fit_null <- sibship_foraging_model(
  colony_count_at_traps, 
  floral_cover_at_traps, 
  trap_coordinates,
  landscape_covariates,
  resistance_model,
  parameter_grid_null,
  verbose=TRUE
)
fit_full <- sibship_foraging_model(
  colony_count_at_traps, 
  floral_cover_at_traps, 
  trap_coordinates,
  landscape_covariates,
  resistance_model,
  parameter_grid_full,
  verbose=TRUE
)
save(fit_null, fit_full, file="midnorth.age_and_road.fitted.RData")

#simulate/refit at maximum likelihood estimates of the parameters
boot_at_mle <- parametric_bootstrap(fit, fit$mle, num_boot=100, verbose=TRUE, random_seed=80)

#simulate/refit at null model
null <- c("theta_age" = mle["theta_age"], "theta_road"=0)
boot_at_null <- parametric_bootstrap(fit, null, num_boot=100, verbose=TRUE, random_seed=80)

save(boot_at_null, boot_at_mle, file="north.stand_age.simulations.RData")

#make some figures to visualize loglik surface, uncertainty in estimates
load("north.stand_age.fitted.RData")
load("north.stand_age.simulations.RData")

dir.create("fig")

plot_1d_likelihood_surface(
  fit, 
  simulations=boot_at_mle, 
  sim_color="dodgerblue"
  ) + 
  xlim(-0.5, 0.5) + 
  ggtitle("Bootstrap") 
ggplot2::ggsave("fig/north.sim_at_mle.png", height=4, width=7, units="in", dpi=300)

plot_1d_likelihood_surface(
  fit, 
  simulations=boot_at_null, 
  sim_color="firebrick"
  ) + 
  xlim(-0.5, 0.5) + 
  ggtitle("Null model") 
ggplot2::ggsave("fig/north.sim_at_null.png", height=4, width=7, units="in", dpi=300)

plot_1d_sampling_distributions(
  fit, 
  parametric_bootstraps=boot_at_mle, 
  null_simulations=boot_at_null, 
  null_color="firebrick", boot_color="dodgerblue"
  ) +
  xlim(-0.5, 0.5) + 
  ggtitle("Null/bootstrap distributions")
ggplot2::ggsave("fig/north.sampling_dist.png", height=4, width=7, units="in", dpi=300)
