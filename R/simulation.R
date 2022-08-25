
simulate_3parameter_model <- function(data, pars, random_seed=NULL)
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
  stopifnot(length(pars) == 3)
  stopifnot(all(names(pars) %in% par_names))
  stopifnot(length(data) == 4)
  stopifnot(all(names(data) %in% data_names))
  pars <- pars[par_names]

  set.seed(random_seed)

  num_colonies <- nrow(data$colony_count_at_traps)
  num_landscape <- length(data$landscape_age)
  num_traps <- ncol(data$colony_count_at_traps)

  # simulate colony locations
  prior_prob_of_nest_location <- 
    prop.table(exp(pars["landscape_age_on_colony_location"] * data$landscape_age))
  colony_locations <- sample(
    1:num_landscape, num_colonies, replace = TRUE, prob = prior_prob_of_nest_location)

  # number of foragers per colony
  foragers_per_colony <- rowSums(data$colony_count_at_traps)

  # simulate forager captures
  floral_cover_at_traps <- outer(data$floral_cover_at_traps, rep(1, num_landscape))
  conditional_capture_rate <-
    prop.table(exp(
      pars["landscape_distance_on_capture_rate"] * data$landscape_distance_to_traps +
      pars["floral_cover_on_capture_rate"] * floral_cover_at_traps), 2)
  colony_count_at_traps <- 
    t(sapply(1:num_colonies, function(i) 
      c(rmultinom(1, size=foragers_per_colony[i], prob=conditional_capture_rate[,colony_locations[i]]))))
  
  #TODO: return intermediates, these might be useful for debugging
  data$colony_count_at_traps <- colony_count_at_traps
  return(data)
}

simulate_2parameter_model <- function(data, pars, random_seed=NULL)
{
  if ("landscape_age_on_colony_location" %in% names(pars))
  {
    pars["landscape_age_on_colony_location"] = 0
  } else {
    pars <- c(pars, "landscape_age_on_colony_location"=0)
  }
  simulate_3parameter_model(data, pars, random_seed)
}

scale_input_data <- function(data){
  data$landscape_distance_to_traps <- 
    (data$landscape_distance_to_traps - mean(data$landscape_distance_to_traps))/sd(data$landscape_distance_to_traps)
  data$floral_cover_at_traps <- 
    (data$floral_cover_at_traps - mean(data$floral_cover_at_traps))/sd(data$floral_cover_at_traps)
  data$landscape_age <- 
    (data$landscape_age - mean(data$landscape_age))/sd(data$landscape_age)
  data
}

