functions {

  real round_nearest(real x, real stepsize){
    return( round(x / stepsize) * stepsize);
  }

}

data {

  ////////////////////////////////
  // general data
  ////////////////////////////////

  int<lower = 0> n_experiments;
  int<lower = 0> n_replicates;
  int<lower = 0> n_total_datapoints;
  vector[n_total_datapoints] detection_limit_log_titer;
  vector<lower =
    min(detection_limit_log_titer)>[n_total_datapoints] observed_log_titers;
  vector<lower = 0>[n_total_datapoints] times;
  int<lower = 0> replicate_id[n_total_datapoints];
  int<lower = 0> experiment_id[n_total_datapoints];

  vector<lower = 0>[n_total_datapoints] n_wells;

  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real intercept_prior_mean;
  real<lower = 0> intercept_prior_sd;

  real decay_rate_prior_mean;
  real<lower = 0> decay_rate_prior_sd;
  
  real<lower = 0> sigma_prior_mean;
  real<lower = 0> sigma_prior_sd;

  real lower_lim_decay_rate;
  
}

transformed data {
  real lld;

  lld = lower_lim_decay_rate;
}


parameters{
  vector<lower = 0>[n_experiments] sigma; // sd about prediction
  vector[n_replicates] intercept[n_experiments];
  vector<lower=lld>[n_experiments] decay_rate;
}

transformed parameters {
  vector[n_total_datapoints] predicted_titer; // predicted titer

  for (i_obs_p in 1:n_total_datapoints){
    int i_exp = experiment_id[i_obs_p];
    int i_repl = replicate_id[i_obs_p];
    
    predicted_titer[i_obs_p] = intercept[i_exp][i_repl] -
      decay_rate[i_exp] * times[i_obs_p];
  }

}

model {

  for (i_obs in 1:n_total_datapoints) {
    real mu = predicted_titer[i_obs];
    real obs_log_titer = observed_log_titers[i_obs];
    real sigma_o = sigma[experiment_id[i_obs]];
    real log_titer_stepsize = 1 / n_wells[i_obs];
    real detection_lim = detection_limit_log_titer[i_obs];

    if (obs_log_titer <= detection_lim) {
      // if this becomes numerically unstable,
      // substitute log(Phi((detection_lim - mu) / sigma))
      target += normal_lcdf(detection_lim |  mu, sigma_o);
    } else {
      target +=
        log_diff_exp(normal_lcdf(obs_log_titer | mu,
                                 sigma_o),
                     normal_lcdf(obs_log_titer - log_titer_stepsize | mu,
                                 sigma_o));
    }
  } // close loop over observations

  // priors
  for (i_exp in 1:n_experiments){
    intercept[i_exp] ~ normal(intercept_prior_mean,
                              intercept_prior_sd);
  }
  
  decay_rate ~ normal(decay_rate_prior_mean,
                 decay_rate_prior_sd);
  sigma ~ normal(sigma_prior_mean,
                 sigma_prior_sd);
}

generated quantities {
  vector[n_total_datapoints] titer_rep_raw;
  vector[n_total_datapoints] titer_rep_censored;

  // posterior predictive checks
  for (i_obs in 1:n_total_datapoints) {
    real dlt = detection_limit_log_titer[i_obs];

    titer_rep_raw[i_obs] = normal_rng(predicted_titer[i_obs],
                                      sigma[experiment_id[i_obs]]);

    titer_rep_censored[i_obs] =
      max({dlt,
            dlt +
            round_nearest(titer_rep_raw[i_obs] - dlt,
                          1.0 / n_wells[i_obs])});
  }
}
