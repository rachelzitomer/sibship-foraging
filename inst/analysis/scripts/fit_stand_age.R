#!/usr/bin/env Rscript
library(optparse)
option_list <- list(
  make_option(c("--region"), type="character", default=NULL, help="One of 'north', 'south', 'midnorth'"),
  make_option(c("--resolution"), type="character", default=NULL, help="One of '30m', '60m', '90m'"),
  make_option(c("--species"), type="character", default=NULL, help="One of 'bomvos', 'bomcal'"),
  make_option(c("--lower_bound"), type="numeric", default=-2.5, help="Lower bound for stand age parameter (scaled by standard deviations)"),
  make_option(c("--upper_bound"), type="numeric", default=2.5, help="Upper bound for stand age parameter (scaled by standard deviations)"),
  make_option(c("--grid_size"), type="integer", default=51, help="Size of stand age parameter grid"),
  make_option(c("--bootstraps"), type="integer", default=0, help="Number of simulations at null model/MLE"),
  make_option(c("--block_size"), type="integer", default=5000, help="Higher might be faster but will use more memory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

region <- opt$region
resolution <- opt$resolution
species <- opt$species
lower_bound <- opt$lower_bound
upper_bound <- opt$upper_bound
grid_size <- opt$grid_size
bootstraps <- opt$bootstraps
block_size <- opt$block_size
prefix <- paste0(region, ".", species, ".", resolution)

#-----------------------------------#
library(sibships)
library(raster)

load(system.file(paste0("data/", prefix, ".RData"), package="sibships"))

landscape_covariates <- raster::stack(list("stand_age"=stand_age))

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

scaling <- 1/sd(values(landscape_covariates[["stand_age"]]))
parameter_grid <- as.matrix(expand.grid(
  "theta"=scaling*seq(lower_bound, upper_bound, length.out=grid_size)
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
save(fit, file=paste0(prefix, ".stand_age.fitted.RData"))

if (bootstraps > 0)
{
  #simulate/refit at maximum likelihood estimates of the parameters
  boot_at_mle <- parametric_bootstrap(
    fit, 
    fit$mle, 
    num_boot=100, 
    verbose=TRUE, 
    random_seed=1
  )
  
  #simulate/refit at null model
  null <- c("theta" = 0)
  boot_at_null <- parametric_bootstrap(
    fit, 
    null, 
    num_boot=100, 
    verbose=TRUE, 
    random_seed=1
  )
  
  save(boot_at_null, boot_at_mle, file=paste0(prefix, ".stand_age.bootstrap.RData"))
}

##make some figures to visualize loglik surface, uncertainty in estimates
##TODO: use prefix in nameing
#
#dir.create("fig")
#
#plot_1d_likelihood_surface(
#  fit, 
#  simulations=boot_at_mle, 
#  sim_color="dodgerblue"
#  ) + 
#  xlim(-0.5, 0.5) + 
#  ggtitle("Bootstrap") 
#ggplot2::ggsave("fig/midnorth.sim_at_mle.png", height=4, width=7, units="in", dpi=300)
#
#plot_1d_likelihood_surface(
#  fit, 
#  simulations=boot_at_null, 
#  sim_color="firebrick"
#  ) + 
#  xlim(-0.5, 0.5) + 
#  ggtitle("Null model") 
#ggplot2::ggsave("fig/midnorth.sim_at_null.png", height=4, width=7, units="in", dpi=300)
#
#plot_1d_sampling_distributions(
#  fit, 
#  parametric_bootstraps=boot_at_mle, 
#  null_simulations=boot_at_null, 
#  null_color="firebrick", boot_color="dodgerblue"
#  ) +
#  xlim(-0.5, 0.5) + 
#  ggtitle("Null/bootstrap distributions")
#ggplot2::ggsave("fig/midnorth.sampling_dist.png", height=4, width=7, units="in", dpi=300)
