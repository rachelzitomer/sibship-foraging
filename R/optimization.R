.collapse_identical_colonies <- function(colony_count_at_traps)
{
  dedup_counts <- unique(colony_count_at_traps)
  A <- table(apply(colony_count_at_traps, 1, paste, collapse="_"))
  B <- apply(dedup_counts, 1, paste, collapse="_")
  weights <- as.numeric(as.matrix(A[match(B, names(A))]))
  return(list(weights=weights, colony_count_at_traps=dedup_counts))
}

fit_3parameter_model <- function(data, start_vals, verbose = TRUE, convergence_tolerance = 1e-4, visitation_always_decreases_with_distance = TRUE)
{
  par_names <- c(
    "floral_cover_on_capture_rate", 
    "landscape_distance_on_capture_rate", 
    "landscape_age_on_colony_location"
  )
  data_names <- c(
    "landscape_distance_to_traps",
    "colony_count_at_traps",
    "floral_cover_at_traps",
    "landscape_age"
  )
  stopifnot(length(start_vals) == 3)
  stopifnot(all(names(start_vals) %in% par_names))
  stopifnot(length(data) == 4)
  stopifnot(all(names(data) %in% data_names))
  upper_bounds <- if (visitation_always_decreases_with_distance) c(Inf, 0, Inf) else c(Inf, Inf, Inf)
  start_vals <- start_vals[par_names]
  ..iter <<- 0
  ..grad <<- c(0, 0, 0)
  #par <- start_vals; browser()#DEBUG
  fit <- nloptr::lbfgs(start_vals, 
    function(par) {
      ..iter <<- ..iter + 1
      ll <- landscape_model_gradient(
        par,
        data$landscape_distance_to_traps, 
        data$colony_count_at_traps, 
        data$floral_cover_at_traps, 
        data$landscape_age
      ) 
      if (verbose) cat("[", ..iter, "] loglik = ", ll[1], "\n", sep="")
      ..grad <<- -ll[-c(1)]
      -ll[1]
    }, gr = function(par) { return(..grad) }, 
    control = list(xtol_rel = convergence_tolerance),
    lower = c(-Inf, -Inf, -Inf),
    upper = upper_bounds
  )
  names(fit$par) <- par_names
  fit
}

fit_2parameter_model_deprecated <- function(data, start_vals, verbose = TRUE, convergence_tolerance = 1e-4, visitation_always_decreases_with_distance = TRUE)
{
  par_names <- c(
    "floral_cover_on_capture_rate", 
    "landscape_distance_on_capture_rate"
  )
  data_names <- c(
    "landscape_distance_to_traps",
    "colony_count_at_traps",
    "floral_cover_at_traps",
    "landscape_age"
  )
  stopifnot(length(start_vals) == 2)
  stopifnot(all(names(start_vals) %in% par_names))
  stopifnot(length(data) == 4)
  stopifnot(all(names(data) %in% data_names))
  upper_bounds <- if (visitation_always_decreases_with_distance) c(Inf, 0) else c(Inf, Inf)
  start_vals <- start_vals[par_names]
  ..iter <<- 0
  ..grad <<- c(0, 0)
  fit <- nloptr::lbfgs(
    start_vals, 
    function(par) {
      ..iter <<- ..iter + 1
      ll <- landscape_model_gradient(
        c(par, 0),
        data$landscape_distance_to_traps, 
        data$colony_count_at_traps, 
        data$floral_cover_at_traps, 
        data$landscape_age
      ) 
      if (verbose) cat("[", ..iter, "] loglik = ", ll[1], "\n", sep="")
      ..grad <<- -ll[-c(1,4)]
      -ll[1]
    }, gr = function(par) { return(..grad) }, 
    control = list(xtol_rel = convergence_tolerance),
    lower = c(-Inf, -Inf),
    upper = upper_bounds
  )
  names(fit$par) <- par_names
  fit
}

fit_2parameter_model <- function(data, start_vals, verbose = TRUE, convergence_tolerance = 1e-4, visitation_always_decreases_with_distance = TRUE)
{
  par_names <- c(
    "floral_cover_on_capture_rate", 
    "landscape_distance_on_capture_rate"
  )
  data_names <- c(
    "landscape_distance_to_traps",
    "colony_count_at_traps",
    "floral_cover_at_traps",
    "landscape_age"
  )
  stopifnot(length(start_vals) == 2)
  stopifnot(all(names(start_vals) %in% par_names))
  stopifnot(length(data) == 4)
  stopifnot(all(names(data) %in% data_names))
  upper_bounds <- if (visitation_always_decreases_with_distance) c(Inf, 0) else c(Inf, Inf)
  start_vals <- start_vals[par_names]
  dedup <- .collapse_identical_colonies(data$colony_count_at_traps)
  ..iter <<- 0
  ..grad <<- c(0, 0)
  fit <- nloptr::lbfgs(
    start_vals, 
    function(par) {
      ..iter <<- ..iter + 1
      ll <- landscape_model_gradient_weighted(
        c(par, 0),
        data$landscape_distance_to_traps, 
        dedup$colony_count_at_traps, 
        data$floral_cover_at_traps, 
        data$landscape_age,
        dedup$weights
      ) 
      if (verbose) cat("[", ..iter, "] loglik = ", ll[1], "\n", sep="")
      ..grad <<- -ll[-c(1,4)]
      -ll[1]
    }, gr = function(par) { return(..grad) }, 
    control = list(xtol_rel = convergence_tolerance),
    lower = c(-Inf, -Inf),
    upper = upper_bounds
  )
  names(fit$par) <- par_names
  fit
}
