#source(system.file("analysis/singletons/midnorth.singletons.R", package="sibships"))

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

#does it work? evaluate at first point in parameter grid
plot(resistance_model(landscape_covariates, parameter_grid[1,]))

#fit_with_singletons <- sibship_foraging_model(
#  colony_count_at_traps, 
#  floral_cover_at_traps, 
#  trap_coordinates,
#  landscape_covariates,
#  resistance_model,
#  parameter_grid,
#  verbose=TRUE
#)
#
#fit_with_no_singletons <- sibship_foraging_model(
#  colony_count_at_traps[rowSums(colony_count_at_traps) > 1,], 
#  floral_cover_at_traps, 
#  trap_coordinates,
#  landscape_covariates,
#  resistance_model,
#  parameter_grid,
#  verbose=TRUE
#)

fit_with_only_singletons <- sibship_foraging_model(
  colony_count_at_traps[rowSums(colony_count_at_traps) == 1,], 
  floral_cover_at_traps, 
  trap_coordinates,
  landscape_covariates,
  resistance_model,
  parameter_grid,
  verbose=TRUE
)

load("midnorth.singletons.RData")
save(fit_with_singletons, fit_with_no_singletons, fit_with_only_singletons, file="midnorth.singletons.RData")

#TODO fig
dir.create("fig")
