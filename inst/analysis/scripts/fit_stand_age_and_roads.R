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
  make_option(c("--block_size"), type="integer", default=25000, help="Higher might be faster but will use more memory")
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

#in practice, we needed to run this on multicore processing, which looked like this:

#wrap bootstrap function so it doesn't require inputs
wrapped_boot<-function(num_boot = 1, random_seed = sample.int(10e6,1)){
  parametric_bootstrap(
    fit, 
    fit$mle, #pars
    num_boot = num_boot,#changed to run one at a time 
    cells_per_block = 25000,
    verbose = TRUE,
    random_seed = random_seed,#different seed for every boot, but actual seeds will be saved in output DF if I need to replicate
    visitation_always_decreases_with_distance = TRUE
  )
}

#run bootstraps
system.time({  
  boot.list<-parallel::mclapply(rep(1,opt$bootstraps), wrapped_boot, mc.cores = 100, mc.preschedule = FALSE)
})

#extract MLEs and loglik_surfaces into dataframes to write out

boot.list.t<-boot.list %>% 
  transpose

boot.df.mles<-as.data.frame(do.call(rbind,boot.list.t$MLEs))

boot.df.loglik.surfaces<-as.data.frame(do.call(rbind, boot.list.t$loglik_surfaces))
