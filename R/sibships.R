.scale_data <- function(data)
{
  var_names <- c(
    "landscape_distance_to_traps", 
    "floral_cover_at_traps", 
    "landscape_age"
  )
  for (var in var_names)
  {
    do_scale <- length(unique(range(data[[var]][]))) == 2
    data[[var]][] <- scale(data[[var]][], center=TRUE, scale=do_scale)
  }
  data
}


.solve_quadratic <- function(X, y, window_size=3, force_interior=TRUE)
{
  # fit a multivariable quadratic model and return maximum
  stopifnot(nrow(X) == length(y))
  X <- as.matrix(X)
  Xp <- poly(X, 2, raw=TRUE)
  argmax <- which.max(y)
  best_X <- X[argmax,,drop=FALSE]
  win_X <- list()
  for(i in colnames(X))
  {
    ux <- unique(X[,i])
    uxo <- order(abs(ux - best_X[,i]))[1:window_size]
    win_X[[i]] <- ux[uxo]
  }
  window <- as.matrix(expand.grid(win_X))
  window <- which(duplicated(rbind(window, X))[-c(1:nrow(window))]) #oof
  y <- y[window]
  Xp <- cbind(1, Xp[window,])
  stopifnot(ncol(Xp) <= length(window))
  colnames(Xp) <- paste0("x", 1:ncol(Xp))
  b <- lm(y ~ 0 + Xp)
  obj <- function(x) {
    Xn <- cbind(1, poly(matrix(x, 1, length(x)), 2, raw=TRUE))
    -Xn %*% coef(b)
  }
  sol <- optim(best_X, fn=obj, method="BFGS")
  ll_max <- -sol$value
  mle <- sol$par
  X_range <- apply(X, 2, range)
  interior <- sapply(1:length(mle),
    function(i) mle[i] >= min(X[,i]) & mle <= max(X[,i])
  )
  out <- matrix(c(mle, ll_max), 1, byrow=TRUE)
  colnames(out) <- c(colnames(X), "loglik")
  #if (force_interior & !all(interior)) { out[] <- NA }
  if (force_interior) { out[!interior] <- out[!interior]*Inf }
  return(out)
}


sibship_foraging_model <- function(
  colony_count_at_traps, 
  floral_cover_at_traps, 
  trap_coordinates, 
  landscape_rasters, 
  resistance_model, 
  parameter_grid, 
  cells_per_block=25000,
  convergence_tolerance=1e-4,
  window_size=3,
  verbose=FALSE)
{
  stopifnot(class(trap_coordinates) == "SpatialPoints")
  stopifnot(class(landscape_rasters) == "RasterStack")
  stopifnot(raster::compareCRS(landscape_covariates, trap_coordinates))
  stopifnot(length(floral_cover_at_traps) == nrow(trap_coordinates))
  stopifnot(ncol(colony_count_at_traps) == nrow(trap_coordinates))

  surface <- radish::conductance_surface(
    landscape_rasters, 
    trap_coordinates,
    directions = 8 #hardcoded, could be settable
  )

  starting_values <- c(
    "floral_cover_on_capture_rate"=0, 
    "landscape_distance_on_capture_rate"=-0.5
  )

  resistance_grid <- array(NA, c(ncell(landscape_rasters), nrow(parameter_grid)))
  resistance_distance_grid <- array(NA, c(length(trap_coordinates), ncell(landscape_rasters), nrow(parameter_grid)))
  fit_list <- list()
  for(i in 1:nrow(parameter_grid))
  {
    cat("Fitting parameters", parameter_grid[i,], "\n", sep=" ")
    resistance_grid[,i] <- values(resistance_model(landscape_rasters, parameter_grid[i,]))
    resistance_distance_grid[,,i] <- distance_to_focal_raw(
      conductance=1./resistance_grid[,i,drop=FALSE],
      s=surface, 
      cells_per_block=cells_per_block,
      average_conductance=FALSE
    )
    data <- list(
      landscape_distance_to_traps = resistance_distance_grid[,,i],
      landscape_age = rep(0, ncell(landscape_rasters)),
      floral_cover_at_traps = floral_cover_at_traps,
      colony_count_at_traps = colony_count_at_traps
    )
    data <- .scale_data(data)
    fit_list[[i]] <- fit_2parameter_model(
      data, 
      starting_values,
      convergence_tolerance=convergence_tolerance,
      verbose=verbose
    )
  }
  fit_ll <- Reduce(c, lapply(fit_list, function(x) -x$value))
  fit_par <- Reduce(rbind, lapply(fit_list, function(x) x$par))
  fitted_models <- cbind(
    parameter_grid, 
    fit_par,
    loglik=fit_ll
  )

  # get better estimate of MLE by quadratic interpolation
  mle <- .solve_quadratic(parameter_grid, fit_ll, window_size=window_size)
  loglik <- mle[1, "loglik"]
  mle <- mle[1, colnames(parameter_grid)]

  #change so that returns -Inf, Inf -- then any(is.infinite(mle)), drop warning
  #if (any(is.na(mle)))
  if (any(is.infinite(mle)))
  {
    #warnings("MLE is outside bounds of parameter grid")
    fit_mle <- resistance <- resistance_distance <- 
      nuisance_parameters <- colony_locations <- NULL
  } else {
    # refit at MLE
    cat("Fitting at MLE\n")
    resistance <- values(resistance_model(landscape_rasters, mle))
    resistance_distance <- distance_to_focal_raw(
      conductance=1./resistance,
      s=surface, 
      cells_per_block=cells_per_block,
      average_conductance=FALSE
    )
    data$landscape_distance_to_traps <- resistance_distance
    data <- .scale_data(data)
    fit_mle <- fit_2parameter_model(
      data, 
      starting_values,
      convergence_tolerance=convergence_tolerance,
      verbose=verbose
    )
    nuisance_parameters <- fit_mle$par
    attr(loglik, "error") <- (-fit_mle$value) - loglik

    # get fitted values at MLE
    colony_locations <- landscape_model_fitted(
      c(nuisance_parameters, 0),
      data$landscape_distance_to_traps, 
      data$colony_count_at_traps, 
      data$floral_cover_at_traps, 
      data$landscape_age
    ) 
  }

  output <- list(
    optimizer=fit_list, 
    parameter_grid=parameter_grid, 
    resistance_grid=resistance_grid,
    resistance_distance_grid=resistance_distance_grid, 
    data=data, 
    landscape_rasters=landscape_rasters,
    surface=surface,
    fitted_models=fitted_models,
    mle=mle,
    loglik=loglik,
    optimizer_mle=fit_mle,
    resistance_distance=resistance_distance,
    resistance=resistance,
    nuisance_parameters=nuisance_parameters,
    colony_locations=colony_locations,
    NULL
  )
  class(output) <- "sibship_foraging_model"
  output
} 


parametric_bootstrap <- function(
  fitted_model, 
  pars, 
  num_boot=10, 
  cells_per_block=25000, 
  window_size=3,
  convergence_tolerance=1e-4,
  random_seed=NULL,
  verbose=FALSE)
{
  stopifnot(class(fitted_model) == "sibship_foraging_model")
  stopifnot(length(pars) == ncol(fitted_model$parameter_grid))
  stopifnot(all(names(pars) %in% colnames(fitted_model$parameter_grid)))

  set.seed(random_seed)

  # Pull out inputs from fitted model object
  resistance_distance_grid <- fitted_model$resistance_distance_grid
  data <- fitted_model$data
  parameter_grid <- fitted_model$parameter_grid
  landscape_rasters <- fitted_model$landscape_rasters
  surface <- fitted_model$surface
  true_parameters <- outer(rep(1, nrow(parameter_grid)), pars)
  stopifnot(all(colnames(true_parameters) == colnames(parameter_grid)))
  colnames(true_parameters) <- paste0("true_", colnames(true_parameters))
  starting_values <- c(
    "floral_cover_on_capture_rate" = 0, 
    "landscape_distance_on_capture_rate" = -0.5
  )

  # Fit model at provided parameter values
  cat("Fitting at provided parameter values\n")
  resistance <- 
    values(resistance_model(landscape_rasters, pars))
  resistance_distance <- distance_to_focal_raw(
    conductance=1./resistance,
    s=surface, 
    cells_per_block=cells_per_block,
    average_conductance=FALSE
  )
  data$landscape_distance_to_traps <- resistance_distance
  data <- .scale_data(data)
  fit <- fit_2parameter_model(
    data, 
    starting_values,
    convergence_tolerance=convergence_tolerance
  )
  nuisance_parameters <- fit$par
  stopifnot(all(names(nuisance_parameters) %in% names(starting_values)))
  true_nuisance_pars <- outer(rep(1, nrow(parameter_grid)), nuisance_parameters)
  colnames(true_nuisance_pars) <- paste0("true_", colnames(true_nuisance_pars))

  # Simulate data from fitted model, then refit to simulations
  seed_sequence <- sample.int(10e6, num_boot)
  fitted_sims <- data.frame()
  mles <- data.frame()
  for(seed in seed_sequence){

    cat("Simulating from bootstrap seed ", seed, " and refitting\n")
    sim_data <- simulate_2parameter_model(
      data, 
      nuisance_parameters, 
      random_seed=seed
    )

    # Fit simulated data across parameter grid
    refit_list <- lapply(1:nrow(parameter_grid), function(i) {
      sim_data$landscape_distance_to_traps <- resistance_distance_grid[,,i]
      sim_data <- .scale_data(sim_data)
      fit_2parameter_model(
        sim_data, 
        starting_values, 
        convergence_tolerance=convergence_tolerance
      )
    })

    # Loglikelihood, nuisance parameters across grid
    sim_loglik <- Reduce(c, lapply(refit_list, function(x) -x$value))
    sim_nuisance_pars <- Reduce(rbind, lapply(refit_list, function(x) x$par))
    colnames(sim_nuisance_pars) <- paste0("sim_", colnames(sim_nuisance_pars))

    # Get better estimate of MLE by quadratic interpolation, refit
    mle <- .solve_quadratic(parameter_grid, sim_loglik, window_size=window_size)
    #if (any(is.na(mle))) 
    #  warning("MLE was outside bounds of parameter grid, setting to NA")

    # Refit simulated data at "true value" (e.g. used to generate simulation)
    sim_data$landscape_distance_to_traps <- resistance_distance
    sim_data <- .scale_data(sim_data)
    fit_truth <- fit_2parameter_model(
      sim_data, 
      starting_values,
      convergence_tolerance=convergence_tolerance,
      verbose=verbose
    )

    # Store MLEs
    tmp <- data.frame(
      seed=seed,
      loglik_at_truth=-fit_truth$value,
      mle
    )
    mles <- rbind(mles, tmp)

    # Store refitted models
    tmp <- cbind(
      seed=seed, 
      true_parameters,
      true_nuisance_pars,
      parameter_grid, 
      sim_nuisance_pars,
      sim_loglik=sim_loglik
    )
    fitted_sims <- rbind(fitted_sims, tmp)
  }

  rownames(fitted_sims) <- NULL
  rownames(mles) <- NULL

  return(list(loglik_surfaces=fitted_sims, MLEs=mles))
} 


nonparametric_bootstrap <- function(
  fitted_model, 
  num_boot=10, 
  window_size=3,
  convergence_tolerance=1e-4,
  random_seed=NULL,
  verbose=FALSE)
{
  stopifnot(class(fitted_model) == "sibship_foraging_model")
  stopifnot(length(pars) == ncol(fitted_model$parameter_grid))
  stopifnot(all(names(pars) %in% colnames(fitted_model$parameter_grid)))

  set.seed(random_seed)

  # Pull out inputs from fitted model object
  resistance_distance_grid <- fitted_model$resistance_distance_grid
  data <- fitted_model$data
  parameter_grid <- fitted_model$parameter_grid
  landscape_rasters <- fitted_model$landscape_rasters
  surface <- fitted_model$surface
  starting_values <- c(
    "floral_cover_on_capture_rate" = 0, 
    "landscape_distance_on_capture_rate" = -0.5
  )

  # Simulate data from fitted model, then refit to simulations
  seed_sequence <- sample.int(10e6, num_boot)
  fitted_sims <- data.frame()
  mles <- data.frame()
  for(seed in seed_sequence){

    cat("Simulating from bootstrap seed ", seed, " and refitting\n")
    sim_data <- data
    bootstrap_sample <- sample(1:nrow(sim_data$colony_count_at_traps), replace=TRUE)
    sim_data$colony_count_at_traps <-
      sim_data$colony_count_at_traps[bootstrap_sample,]

    # Fit simulated data across parameter grid
    refit_list <- lapply(1:nrow(parameter_grid), function(i) {
      sim_data$landscape_distance_to_traps <- resistance_distance_grid[,,i]
      sim_data <- .scale_data(sim_data)
      fit_2parameter_model(
        sim_data, 
        starting_values, 
        convergence_tolerance=convergence_tolerance
      )
    })

    # Loglikelihood, nuisance parameters across grid
    sim_loglik <- Reduce(c, lapply(refit_list, function(x) -x$value))
    sim_nuisance_pars <- Reduce(rbind, lapply(refit_list, function(x) x$par))
    colnames(sim_nuisance_pars) <- paste0("sim_", colnames(sim_nuisance_pars))

    # Get better estimate of MLE by quadratic interpolation, refit
    mle <- .solve_quadratic(parameter_grid, sim_loglik, window_size=window_size)
    if (any(is.na(mle))) warning("MLE was outside bounds of parameter grid, setting to NA")

    # Store MLEs
    tmp <- data.frame(
      seed=seed,
      mle
    )
    mles <- rbind(mles, tmp)

    # Store refitted models
    parameter_grid_rename <- parameter_grid
    colnames(parameter_grid_rename) <-
      paste0("sim_", colnames(parameter_grid_rename))
    tmp <- cbind(
      seed=seed, 
      parameter_grid_rename, 
      sim_nuisance_pars,
      sim_loglik=sim_loglik
    )
    fitted_sims <- rbind(fitted_sims, tmp)
  }

  rownames(fitted_sims) <- NULL
  rownames(mles) <- NULL

  return(list(loglik_surfaces=fitted_sims, MLEs=mles))
} 
