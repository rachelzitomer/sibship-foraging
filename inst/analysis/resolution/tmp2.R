region <- "midnorth"
resolution <- "90m"
species <- "bomvos"
lower_bound <- -2.5
upper_bound <- 2.5
grid_size <- 5
bootstraps <- 10
block_size <- 10000
prefix <- paste0(region, ".", species, ".", resolution)

#-----------------------------------#
library(sibships)
library(raster)

load(system.file(paste0("data/", prefix, ".RData"), package="sibships"))

#TODO:
# it seems like some values of roads raster are NA?
# we replace with 0 here, but this really should be fixed
roads[is.na(roads[])] <- 0

landscape_covariates <- raster::stack(list("stand_age"=stand_age, "roads"=roads))

resistance_model <- function(raster_stack, parameters)
{
  # converts stack of rasters to a single output raster
  # representing a resistance surface
  stopifnot("stand_age" %in% names(raster_stack))
  stopifnot("roads" %in% names(raster_stack))
  stopifnot("stand_age" %in% names(parameters))
  stopifnot("roads" %in% names(parameters))
  stopifnot(length(parameters) == 2)
  resistance <- exp(raster_stack[["stand_age"]] * parameters["stand_age"] + raster_stack[["roads"]] * parameters["roads"])
  return(resistance)
}

browser()

scaling_stand_age <- 1/sd(values(landscape_covariates[["stand_age"]]))
scaling_roads <- 1/sd(values(landscape_covariates[["roads"]]))
parameter_grid <- as.matrix(expand.grid(
  "stand_age"=scaling_stand_age*seq(lower_bound, upper_bound, length.out=grid_size),
  "roads"=scaling_roads*seq(lower_bound, upper_bound, length.out=grid_size)
))

#does model work? evaluate at first point in parameter grid, check for NAs, etc
print(resistance_model(landscape_covariates, parameter_grid[1,]))

fit <- sibship_foraging_model(
  colony_count_at_traps, 
  floral_cover_at_traps, 
  trap_coordinates,
  landscape_covariates,
  resistance_model,
  parameter_grid,
  verbose=TRUE,
  cells_per_block=block_size
)
save(fit, file=paste0(prefix, ".stand_age_and_roads.fitted.RData"))

if (bootstraps > 0)
{
  #simulate/refit at maximum likelihood estimates of the parameters
  boot_at_mle <- parametric_bootstrap(
    fit, 
    fit$mle, 
    num_boot=5, 
    verbose=TRUE, 
    random_seed=1,
    visitation_always_decreases_with_distance=visitation_always_decreases_with_distance
  )
  
  ##simulate/refit at null model
  ##REMOVED: this often ends up on the boundary 
  ##  (e.g. fitted nuisance parameters where there's no effect of landscape on movement)
  ##  in which case there's no MLE for the resistance distance parameters!
  boot_at_null <- NA
  #null <- c("theta" = 0)
  #boot_at_null <- parametric_bootstrap(
  #  fit, 
  #  null, 
  #  num_boot=100, 
  #  verbose=TRUE, 
  #  random_seed=1,
  #  visitation_always_decreases_with_distance=visitation_always_decreases_with_distance
  #)
  
  save(boot_at_null, boot_at_mle, file=paste0(prefix, ".stand_age_and_roads.bootstrap.RData"))
}

