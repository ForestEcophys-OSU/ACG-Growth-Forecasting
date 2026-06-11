data {
  int<lower=1> N_obs; // Total number of observations
  vector[N_obs] y_obs; // Observed DBHs (on original scale)

  int<lower=1> N_trees; // Number of individual trees
  int<lower=1> N_unique_years; // Number of unique years

  array[N_obs] int<lower=1, upper=N_unique_years> year_idx_obs; // Index of year for each observation
  array[N_obs] int<lower=1, upper=N_trees> tree_idx_obs; // Index of tree for each observation

  vector[N_unique_years] unique_years; // Sorted unique years in the data

  // For initial conditions of latent states (log-scale)
  real log_x_ic_mean; // Mean of initial log(DBH)
  real<lower=0> log_x_ic_sd; // SD of initial log(DBH)
}

parameters {
  // Raw latent states for non-centered parameterization (on log scale)
  matrix[N_unique_years, N_trees] log_x_raw;

  // Standard deviations (on original and log scales)
  real<lower=0> sigma_obs; // Observation error (on original DBH scale)
  real<lower=0> sigma_add_log; // Process error (on log(DBH) scale)

  // Growth model parameters
  real<lower=0> R; // Growth rate (constrained positive)
  real<lower=0> K; // Carrying capacity (constrained positive)
}

transformed parameters {
  // True latent states (on log scale, derived from raw and process model)
  matrix[N_unique_years, N_trees] log_x;
  // True latent states (on original scale)
  matrix<lower=0>[N_unique_years, N_trees] x;

  // Initial condition for latent states (log_x[1,r])
  for (r in 1:N_trees) {
    log_x[1, r] = log_x_ic_mean + log_x_ic_sd * log_x_raw[1, r];
  }

  // Process model loop
  for (t_idx in 2:N_unique_years) {
    real dt = unique_years[t_idx] - unique_years[t_idx - 1];
    for (r in 1:N_trees) {
      // Current DBH on original scale
      real current_x = exp(log_x[t_idx - 1, r]);

      // Expected absolute growth increment (logistic model)
      real growth_increment = R * current_x * (1 - current_x / K) * dt;

      // Expected next true DBH (on original scale)
      real expected_next_x_original = current_x + growth_increment;

      // Ensure expected_next_x_original is positive before taking log
      // This guard is crucial for numerical stability
      if (expected_next_x_original <= 0) {
        expected_next_x_original = 0.001; // Small positive value, slightly larger to give more room.
      }
      
      // Expected next true log(DBH)
      real expected_log_next_x = log(expected_next_x_original);

      // Non-centered reparameterization for process error
      log_x[t_idx, r] = expected_log_next_x + sigma_add_log * log_x_raw[t_idx, r];
    }
  }

  // Convert log_x back to original scale (x) for observation model
  for (t in 1:N_unique_years) {
    for (r in 1:N_trees) {
      x[t, r] = exp(log_x[t, r]);
    }
  }
}

model {
  // Priors for standard deviations (using Half-Cauchy for more robustness)
  sigma_obs ~ cauchy(0, 1);     // Half-Cauchy(0, scale=1) is a common robust prior for std dev
  sigma_add_log ~ cauchy(0, 0.1); // Half-Cauchy(0, scale=0.1) for log-scale process error

  // Priors for growth model parameters (more informative priors based on simulation)
  R ~ lognormal(log(0.1), 0.2); // Mean around 0.1, std dev of log(R) = 0.2 (allows R to range from ~0.06 to ~0.15)
  K ~ normal(60, 10) T[0,];    // Mean around 60 cm, allowing values from ~30-90 cm (more constrained than 20)

  // Priors for raw latent states (standard normal)
  to_vector(log_x_raw) ~ normal(0, 1);

  // Observation model loop
  for (i in 1:N_obs) {
    y_obs[i] ~ normal(x[year_idx_obs[i], tree_idx_obs[i]], sigma_obs);
  }
}

generated quantities {
  vector[N_obs] y_pred_mean;

  for (i in 1:N_obs) {
    y_pred_mean[i] = x[year_idx_obs[i], tree_idx_obs[i]];
  }
  matrix[N_unique_years, N_trees] x_mean_output = x;
}

