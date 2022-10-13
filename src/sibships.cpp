#include <RcppArmadillo.h> 

// TODO check gradient with nonzero offsets

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::vec softmax_jacobian_multiply (arma::vec gradient, arma::vec z)
{
  if (gradient.n_elem != z.n_elem) Rcpp::stop("Dimension mismatch");
  arma::vec out (gradient.n_elem, arma::fill::zeros);
  double gr_z = arma::accu(gradient % z);
  for (unsigned i=0; i<gradient.n_elem; ++i)
  {
    //for (unsigned j=0; j<gradient.n_elem; ++j)
    //{
    //  out.at(i) += z.at(i) * (int(i==j) - z.at(j)) * gradient.at(j);
    //}
    // simplifies to: -out.at(i) = z.at(i) * sum(z * gradient) - z.at(i) * gradient.at(i)
    out.at(i) = -z.at(i) * gr_z + z.at(i) * gradient.at(i);
  }
  return out;
}

struct landscape_model
{
  // evaluate loglikelihood of basic captures x landscape model
  
  const unsigned num_landscape;
  const unsigned num_colonies;
  const unsigned num_traps;
  const arma::mat landscape_distance_to_traps;
  const arma::mat colony_counts_at_traps;
  const arma::vec floral_cover_at_traps;
  const arma::vec landscape_age;
  const arma::vec weights;
  const arma::vec offsets;

  arma::mat conditional_capture_rate;
  arma::vec prior_prob_of_nest_location;

  // DEPRECATE ME
  landscape_model (arma::mat landscape_distance_to_traps,
                   arma::mat colony_counts_at_traps,
                   arma::vec floral_cover_at_traps,
                   arma::vec landscape_age)
    : landscape_distance_to_traps (landscape_distance_to_traps)
    , colony_counts_at_traps (colony_counts_at_traps)
    , floral_cover_at_traps (floral_cover_at_traps)
    , landscape_age (landscape_age)
    , num_traps (floral_cover_at_traps.n_elem)
    , num_colonies (colony_counts_at_traps.n_rows)
    , num_landscape (landscape_age.n_elem)
    , weights (arma::ones(num_colonies))
    , offsets (arma::zeros(num_traps))
    , conditional_capture_rate (num_traps, num_landscape)
    , prior_prob_of_nest_location (num_landscape)
  {
    if(!(landscape_distance_to_traps.n_rows == num_traps && landscape_distance_to_traps.n_cols == num_landscape)) Rcpp::stop("err dim(landscape_distance_to_traps)");
    if(!(colony_counts_at_traps.n_rows == num_colonies && colony_counts_at_traps.n_cols == num_traps)) Rcpp::stop("err dim(colony_counts_at_traps)");
  }

  // DEPRECATE ME
  landscape_model (arma::mat landscape_distance_to_traps,
                   arma::mat colony_counts_at_traps,
                   arma::vec floral_cover_at_traps,
                   arma::vec landscape_age,
                   arma::vec weights)
    : landscape_distance_to_traps (landscape_distance_to_traps)
    , colony_counts_at_traps (colony_counts_at_traps)
    , floral_cover_at_traps (floral_cover_at_traps)
    , landscape_age (landscape_age)
    , num_traps (floral_cover_at_traps.n_elem)
    , num_colonies (colony_counts_at_traps.n_rows)
    , num_landscape (landscape_age.n_elem)
    , weights (weights)
    , offsets (arma::zeros(num_traps))
    , conditional_capture_rate (num_traps, num_landscape)
    , prior_prob_of_nest_location (num_landscape)
  {
    if(!(landscape_distance_to_traps.n_rows == num_traps && landscape_distance_to_traps.n_cols == num_landscape)) Rcpp::stop("err dim(landscape_distance_to_traps)");
    if(!(colony_counts_at_traps.n_rows == num_colonies && colony_counts_at_traps.n_cols == num_traps)) Rcpp::stop("err dim(colony_counts_at_traps)");
    if(!(weights.n_elem == num_colonies)) Rcpp::stop("err dim(weights)");
  }

  landscape_model (arma::mat landscape_distance_to_traps,
                   arma::mat colony_counts_at_traps,
                   arma::vec floral_cover_at_traps,
                   arma::vec landscape_age,
                   arma::vec weights,
                   arma::vec offsets)
    : landscape_distance_to_traps (landscape_distance_to_traps)
    , colony_counts_at_traps (colony_counts_at_traps)
    , floral_cover_at_traps (floral_cover_at_traps)
    , landscape_age (landscape_age)
    , num_traps (floral_cover_at_traps.n_elem)
    , num_colonies (colony_counts_at_traps.n_rows)
    , num_landscape (landscape_age.n_elem)
    , weights (weights)
    , offsets (offsets)
    , conditional_capture_rate (num_traps, num_landscape)
    , prior_prob_of_nest_location (num_landscape)
  {
    if(!(landscape_distance_to_traps.n_rows == num_traps && landscape_distance_to_traps.n_cols == num_landscape)) Rcpp::stop("err dim(landscape_distance_to_traps)");
    if(!(colony_counts_at_traps.n_rows == num_colonies && colony_counts_at_traps.n_cols == num_traps)) Rcpp::stop("err dim(colony_counts_at_traps)");
    if(!(weights.n_elem == num_colonies)) Rcpp::stop("err dim(weights)");
    if(!(offsets.n_elem == num_traps)) Rcpp::stop("err dim(offsets)");
  }

  double marginal_likelihood (const double floral_cover_on_capture_rate, const double landscape_distance_on_capture_rate, const double landscape_age_on_colony_location)
  {
    // conditional capture rates
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      for (unsigned trap=0; trap<num_traps; ++trap)
      {
        conditional_capture_rate(trap, cell) = offsets(trap) +
          floral_cover_on_capture_rate * floral_cover_at_traps(trap) +
              landscape_distance_on_capture_rate * landscape_distance_to_traps(trap, cell);
      }
      conditional_capture_rate.col(cell) -= conditional_capture_rate.col(cell).max();
      conditional_capture_rate.col(cell)  = arma::exp(conditional_capture_rate.col(cell));
      conditional_capture_rate.col(cell) /= arma::accu(conditional_capture_rate.col(cell));
    }

    // prior on colony location
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      prior_prob_of_nest_location(cell) = landscape_age_on_colony_location * landscape_age(cell);
    }
    prior_prob_of_nest_location -= prior_prob_of_nest_location.max();
    prior_prob_of_nest_location  = arma::exp(prior_prob_of_nest_location);
    prior_prob_of_nest_location /= arma::accu(prior_prob_of_nest_location);

    // constants
    arma::vec constants (num_colonies, arma::fill::zeros);
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      for (unsigned trap=0; trap<num_traps; ++trap) // multinomial loglikelihood
      {
        constants.at(colony) -= R::lgammafn(colony_counts_at_traps.at(colony,trap) + 1.);
      }
      constants.at(colony) += R::lgammafn(arma::accu(colony_counts_at_traps.row(colony)) + 1.);
    }

    // marginalize
    double target = 0.;
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      double marginal_lik = 0.;

      for (unsigned cell=0; cell<num_landscape; ++cell)
      {
        double conditional_lik = log(prior_prob_of_nest_location.at(cell)); 

        for (unsigned trap=0; trap<num_traps; ++trap) // multinomial loglikelihood
        {
          conditional_lik += log(conditional_capture_rate.at(trap,cell)) * colony_counts_at_traps.at(colony,trap);
        }

        marginal_lik += exp(conditional_lik + constants.at(colony));
      }

      target += weights.at(colony) * log(marginal_lik);
    }

    return target;
  }

  arma::vec gradient (const double floral_cover_on_capture_rate, const double landscape_distance_on_capture_rate, const double landscape_age_on_colony_location)
  {
    // conditional capture rates
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      for (unsigned trap=0; trap<num_traps; ++trap)
      {
        conditional_capture_rate(trap, cell) = offsets(trap) +
          floral_cover_on_capture_rate * floral_cover_at_traps(trap) +
              landscape_distance_on_capture_rate * landscape_distance_to_traps(trap, cell);
      }
      conditional_capture_rate.col(cell) -= conditional_capture_rate.col(cell).max();
      conditional_capture_rate.col(cell)  = arma::exp(conditional_capture_rate.col(cell));
      conditional_capture_rate.col(cell) /= arma::accu(conditional_capture_rate.col(cell));
    }

    // prior on colony location
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      prior_prob_of_nest_location(cell) = landscape_age_on_colony_location * landscape_age(cell);
    }
    prior_prob_of_nest_location -= prior_prob_of_nest_location.max();
    prior_prob_of_nest_location  = arma::exp(prior_prob_of_nest_location);
    prior_prob_of_nest_location /= arma::accu(prior_prob_of_nest_location);

    // constants
    arma::vec constants (num_colonies, arma::fill::zeros);
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      for (unsigned trap=0; trap<num_traps; ++trap) // multinomial loglikelihood
      {
        constants.at(colony) -= R::lgammafn(colony_counts_at_traps.at(colony,trap) + 1.);
      }
      constants.at(colony) += R::lgammafn(arma::accu(colony_counts_at_traps.row(colony)) + 1.);
    }

    // marginalize
    arma::vec marginal_lik = arma::zeros<arma::vec>(num_colonies);
    arma::mat conditional_lik = arma::zeros<arma::mat>(num_landscape,num_colonies);
    double target = 0.;
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      marginal_lik(colony) = 0.;

      for (unsigned cell=0; cell<num_landscape; ++cell)
      {
        conditional_lik(cell,colony) = log(prior_prob_of_nest_location.at(cell)); 

        for (unsigned trap=0; trap<num_traps; ++trap) // multinomial loglikelihood
        {
          conditional_lik(cell,colony) += log(conditional_capture_rate.at(trap,cell)) * colony_counts_at_traps.at(colony,trap);
        }

        marginal_lik(colony) += exp(conditional_lik(cell,colony) + constants.at(colony));
      }

      target += weights(colony) * log(marginal_lik(colony));
    }

    // marginal likelihood, reverse diff
    arma::vec d_prior_prob_of_nest_location = arma::zeros<arma::vec>(num_landscape);
    arma::mat d_conditional_capture_rate = arma::zeros<arma::mat>(num_traps,num_landscape);
   
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      double d_marginal_lik = weights(colony)/marginal_lik(colony);
      for (unsigned cell=0; cell<num_landscape; ++cell)
      {
        double d_conditional_lik = d_marginal_lik * exp(conditional_lik(cell,colony) + constants.at(colony));
        for (unsigned trap=0; trap<num_traps; ++trap)
        {
          d_conditional_capture_rate(trap,cell) += d_conditional_lik * colony_counts_at_traps(colony,trap) / conditional_capture_rate(trap,cell);
        }
        d_prior_prob_of_nest_location(cell) += d_conditional_lik / prior_prob_of_nest_location(cell);
      }
    }
   
    // prior on colony location, reverse diff
    d_prior_prob_of_nest_location = softmax_jacobian_multiply(d_prior_prob_of_nest_location,  prior_prob_of_nest_location);
    double d_landscape_age_on_colony_location = arma::accu(landscape_age % d_prior_prob_of_nest_location);

    // conditional capture rates, reverse diff
    double d_floral_cover_on_capture_rate = 0;
    double d_landscape_distance_on_capture_rate = 0;
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      d_conditional_capture_rate.col(cell) = softmax_jacobian_multiply(d_conditional_capture_rate.col(cell), conditional_capture_rate.col(cell));
      for (unsigned trap=0; trap<num_traps; ++trap)
      {
        d_floral_cover_on_capture_rate += d_conditional_capture_rate(trap, cell) * floral_cover_at_traps(trap);
        d_landscape_distance_on_capture_rate += d_conditional_capture_rate(trap, cell) * landscape_distance_to_traps(trap, cell);
      }
    }
    
    // [loglik, gradient]
    return arma::vec({target, d_floral_cover_on_capture_rate, d_landscape_distance_on_capture_rate, d_landscape_age_on_colony_location});
  }

  arma::mat conditional_distribution (const double floral_cover_on_capture_rate, const double landscape_distance_on_capture_rate, const double landscape_age_on_colony_location)
  {
    arma::mat conditional_posterior (num_landscape, num_colonies);

    // conditional capture rates
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      for (unsigned trap=0; trap<num_traps; ++trap)
      {
        conditional_capture_rate(trap, cell) = offsets(trap) +
          floral_cover_on_capture_rate * floral_cover_at_traps(trap) +
              landscape_distance_on_capture_rate * landscape_distance_to_traps(trap, cell);
      }
      conditional_capture_rate.col(cell) -= conditional_capture_rate.col(cell).max();
      conditional_capture_rate.col(cell)  = arma::exp(conditional_capture_rate.col(cell));
      conditional_capture_rate.col(cell) /= arma::accu(conditional_capture_rate.col(cell));
    }

    // prior on colony location
    for (unsigned cell=0; cell<num_landscape; ++cell)
    {
      prior_prob_of_nest_location(cell) = landscape_age_on_colony_location * landscape_age(cell);
    }
    prior_prob_of_nest_location -= prior_prob_of_nest_location.max();
    prior_prob_of_nest_location  = arma::exp(prior_prob_of_nest_location);
    prior_prob_of_nest_location /= arma::accu(prior_prob_of_nest_location);

    // constants
    arma::vec constants (num_colonies, arma::fill::zeros);
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      for (unsigned trap=0; trap<num_traps; ++trap) // multinomial loglikelihood
      {
        constants.at(colony) -= R::lgammafn(colony_counts_at_traps.at(colony,trap) + 1.);
      }
      constants.at(colony) += R::lgammafn(arma::accu(colony_counts_at_traps.row(colony)) + 1.);
    }

    // marginalize
    double target = 0.;
    for (unsigned colony=0; colony<num_colonies; ++colony)
    {
      for (unsigned cell=0; cell<num_landscape; ++cell)
      {
        double conditional_lik = log(prior_prob_of_nest_location.at(cell)); 

        for (unsigned trap=0; trap<num_traps; ++trap) // multinomial loglikelihood
        {
          conditional_lik += log(conditional_capture_rate.at(trap,cell)) * colony_counts_at_traps.at(colony,trap);
        }

        conditional_posterior.at(cell, colony) = conditional_lik;
      }

      conditional_posterior.col(colony) -=  conditional_posterior.col(colony).max();
      conditional_posterior.col(colony) = arma::exp(conditional_posterior.col(colony));
      conditional_posterior.col(colony) /= arma::accu(conditional_posterior.col(colony));
    }

    return conditional_posterior;
  }
  
  //arma::vec posterior_predictive (const double floral_cover_on_capture_rate, const double landscape_distance_on_capture_rate, const double landscape_age_on_colony_location)
  //{
  //  //combine these?
  //  arma::mat conditional_posterior = ...;//TODO
  //  arma::mat conditional_capture_rate = ...;//TODO
  //  // loop over colonies, calculate average pairwise distance
  //  arma::vec mean_pairwise_distance(num_colonies);
  //  for(colony=0; colony<?; ++colony)
  //  {
  //    for(unsigned i=0; i<num_traps; ++i)
  //    {
  //      for(unsigned j=0; j<num_traps; ++j)
  //      {
  //        mean_pairwise_distance(colony) += 
  //          conditional_posterior.at(cell, colony) *
  //          conditional_capture_rate.at(i, cell) *
  //          conditional_capture_rate.at(j, cell) *
  //          trap_distances.at(i, j);
  //      }
  //    }
  //  }
  //  return mean_pairwise_distance;
  //}
};

// [[Rcpp::export]]
double landscape_model_loglik (
    arma::vec par, 
    arma::mat landscape_distance_to_traps,
    arma::mat colony_counts_at_traps,
    arma::vec floral_cover_at_traps,
    arma::vec landscape_age
    )
{
  if(par.n_elem != 3) Rcpp::stop("err dim(par)");
  landscape_model model (landscape_distance_to_traps, colony_counts_at_traps, floral_cover_at_traps, landscape_age);
  return model.marginal_likelihood(par[0], par[1], par[2]);
}

// [[Rcpp::export]]
arma::vec landscape_model_gradient (
    arma::vec par, 
    arma::mat landscape_distance_to_traps,
    arma::mat colony_counts_at_traps,
    arma::vec floral_cover_at_traps,
    arma::vec landscape_age
    )
{
  if(par.n_elem != 3) Rcpp::stop("err dim(par)");
  landscape_model model (landscape_distance_to_traps, colony_counts_at_traps, floral_cover_at_traps, landscape_age);
  return model.gradient(par[0], par[1], par[2]);
}

// [[Rcpp::export]]
Rcpp::List landscape_model_fitted (
    arma::vec par, 
    arma::mat landscape_distance_to_traps,
    arma::mat colony_counts_at_traps,
    arma::vec floral_cover_at_traps,
    arma::vec landscape_age
    )
{
  if(par.n_elem != 3) Rcpp::stop("err dim(par)");
  landscape_model model (landscape_distance_to_traps, colony_counts_at_traps, floral_cover_at_traps, landscape_age);
  arma::mat conditional_posterior = model.conditional_distribution(par[0], par[1], par[2]);
  return Rcpp::List::create(
      Rcpp::_["conditional_capture_rate"] = model.conditional_capture_rate,
      Rcpp::_["conditional_posterior"] = conditional_posterior
  );
}

// [[Rcpp::export]]
arma::vec landscape_model_gradient_weighted (
    arma::vec par, 
    arma::mat landscape_distance_to_traps,
    arma::mat colony_counts_at_traps,
    arma::vec floral_cover_at_traps,
    arma::vec landscape_age,
    arma::vec weights
    )
{
  if(par.n_elem != 3) Rcpp::stop("err dim(par)");
  landscape_model model (landscape_distance_to_traps, colony_counts_at_traps, floral_cover_at_traps, landscape_age, weights);
  return model.gradient(par[0], par[1], par[2]);
}

//RCPP_EXPOSED_CLASS_NODECL(landscape_model)
//
//RCPP_MODULE(landscape_model) {
//  using namespace Rcpp;
//  class_<landscape_model>("landscape_model")
//    .constructor<arma::mat, arma::mat, arma::vec, arma::vec>()
//    .method("marginal_likelihood", &landscape_model::marginal_likelihood)
//    ;
//}
