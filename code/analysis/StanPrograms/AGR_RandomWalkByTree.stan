// AGR_RandomWalkByTree_Clean.stan
// A hierarchical random walk state-space model for DBH of multiple trees
// Uses non-centered parameterization for better sampling efficiency

data {
  int<lower=1> num_obs;          // Total number of observations across all trees
  int<lower=1> num_trees;        // Number of individual trees

  array[num_obs] int<lower=1, upper=num_trees> tree_id; // Tree ID for each observation
  vector[num_obs] years_all;     // Years of observation for all data points
  vector[num_obs] observed_dbh_all; // Observed DBH values for all data points

  // Initial true DBH for each tree (for their respective priors)
  vector[num_trees] initial_dbh_mean_prior_per_tree;
  real<lower=0> initial_dbh_sd_prior; // Assuming same SD for all initial tree priors

  // To map observations to time steps within each tree's sequence
  array[num_trees] int<lower=1> N_years_per_tree; // Number of years observed for each tree
  int<lower=1> max_N_years_per_tree;             // Max years observed for any tree

  array[num_trees] vector[max_N_years_per_tree] years_per_tree; // Years for each tree
  array[num_trees] vector[max_N_years_per_tree] observed_dbh_per_tree; // DBH for each tree
}

parameters {
  // Population-level parameters for the hierarchical drift
  real pop_AGR_drift_mean;          // Mean of the population-level AGR_drift
  real<lower=0> pop_AGR_drift_sd;   // Standard deviation of individual AGR_drifts around the population mean

  // Raw individual-level AGR_drift parameters (non-centered)
  vector[num_trees] AGR_drift_raw; // Standard normal deviates for individual AGRs

  // Raw initial true DBH for each tree (non-centered)
  vector[num_trees] initial_dbh_raw; // Standard normal deviates for initial DBH of each tree

  // Raw latent states for non-centered parameterization of the process model (for ALL time points/trees, excluding initial)
  // The size is num_obs - num_trees because the first state for each tree is handled by initial_dbh_raw
  vector[num_obs - num_trees] dbh_process_error_raw; // Standard normal deviates for process error

  // Standard deviations for error components (shared across all trees)
  real<lower=0> sigma_measurement_error; // Measurement error standard deviation
  real<lower=0> sigma_process_error;    // Process error standard deviation
}

transformed parameters {
  // Individual-level AGR_drift (centered from raw)
  vector[num_trees] AGR_drift_per_tree;
  for (j in 1:num_trees) {
    AGR_drift_per_tree[j] = pop_AGR_drift_mean + pop_AGR_drift_sd * AGR_drift_raw[j];
  }
}

model {
  // Local variables for true DBH states and indexing
  vector<lower=0>[num_obs] dbh_true_all; // Combined vector for all true DBH states, local to model block
  int current_obs_idx_global = 1;      // Global index for observed_dbh_all and dbh_true_all
  int current_process_error_idx = 1;   // Index for dbh_process_error_raw

  // Priors for population-level parameters
  pop_AGR_drift_mean ~ normal(0.5, 0.2); // Expect population mean AGR around 0.5 cm/year
  pop_AGR_drift_sd   ~ cauchy(0, 0.1);  // How much individuals vary around the mean AGR

  // Priors for standard deviations (using Half-Cauchy for robustness)
  sigma_measurement_error ~ cauchy(0, 0.1); // Assuming precise measurements (e.g., 1mm error)
  sigma_process_error   ~ cauchy(0, 0.05); // Expecting true year-to-year variation to be small

  // Priors for raw individual-level AGRs (standard normal)
  AGR_drift_raw ~ normal(0, 1);

  // Priors for raw initial DBH values (standard normal)
  initial_dbh_raw ~ normal(0, 1);

  // Priors for the raw process errors (standard normal)
  dbh_process_error_raw ~ normal(0, 1);

  // Loop through each tree to define its true DBH trajectory and likelihood
  for (j in 1:num_trees) {
    // Initial true DBH state for tree j (Non-centered)
    dbh_true_all[current_obs_idx_global] = initial_dbh_mean_prior_per_tree[j] +
                                            initial_dbh_sd_prior * initial_dbh_raw[j];

    // Explicitly enforce lower bound after transformation
    if (dbh_true_all[current_obs_idx_global] < 0.01) {
      dbh_true_all[current_obs_idx_global] = 0.01;
    }

    // Process Model for tree j: Random Walk with Drift (Non-Centered)
    for (t_idx_local in 2:N_years_per_tree[j]) {
      real dt = years_per_tree[j, t_idx_local] - years_per_tree[j, t_idx_local - 1];
      real expected_dbh_next = dbh_true_all[current_obs_idx_global + t_idx_local - 2] + AGR_drift_per_tree[j] * dt;

      // Non-centered reparameterization for process error
      dbh_true_all[current_obs_idx_global + t_idx_local - 1] = expected_dbh_next +
                                                                 sigma_process_error * dbh_process_error_raw[current_process_error_idx];

      // Explicitly enforce positivity for subsequent states
      if (dbh_true_all[current_obs_idx_global + t_idx_local - 1] < 0.01) {
        dbh_true_all[current_obs_idx_global + t_idx_local - 1] = 0.01;
      }
      current_process_error_idx += 1; // Increment for each process error term
    }
    current_obs_idx_global += N_years_per_tree[j]; // Advance global index for next tree
  }

  // Data Model (Observation Model): Linking observed DBH to true DBH
  observed_dbh_all ~ normal(dbh_true_all, sigma_measurement_error); // Vectorized form
}

generated quantities {
  // Local variables for true DBH states and indexing to recreate dbh_true_all
  // We need to re-calculate dbh_true_all here if we want to extract it,
  // as it was a local variable in the model block.
  // This is a common pattern for state-space models.
  vector<lower=0>[num_obs] predicted_dbh_true_all;
  int current_obs_idx_global_gq = 1;
  int current_process_error_idx_gq = 1;

  for (j in 1:num_trees) {
    predicted_dbh_true_all[current_obs_idx_global_gq] = initial_dbh_mean_prior_per_tree[j] +
                                                        initial_dbh_sd_prior * initial_dbh_raw[j];
    if (predicted_dbh_true_all[current_obs_idx_global_gq] < 0.01) {
      predicted_dbh_true_all[current_obs_idx_global_gq] = 0.01;
    }

    for (t_idx_local in 2:N_years_per_tree[j]) {
      real dt = years_per_tree[j, t_idx_local] - years_per_tree[j, t_idx_local - 1];
      real expected_dbh_next = predicted_dbh_true_all[current_obs_idx_global_gq + t_idx_local - 2] + AGR_drift_per_tree[j] * dt;

      predicted_dbh_true_all[current_obs_idx_global_gq + t_idx_local - 1] = expected_dbh_next +
                                                                             sigma_process_error * dbh_process_error_raw[current_process_error_idx_gq];
      if (predicted_dbh_true_all[current_obs_idx_global_gq + t_idx_local - 1] < 0.01) {
        predicted_dbh_true_all[current_obs_idx_global_gq + t_idx_local - 1] = 0.01;
      }
      current_process_error_idx_gq += 1;
    }
    current_obs_idx_global_gq += N_years_per_tree[j];
  }


  // Calculate annual absolute growth rates (AGRs) for each tree from estimated true DBH
  int total_growth_intervals = 0;
  for (j in 1:num_trees) {
    if (N_years_per_tree[j] > 1) {
      total_growth_intervals += (N_years_per_tree[j] - 1);
    }
  }

  vector[total_growth_intervals] all_annual_agr_true;
  int current_idx_agr = 1;
  int current_obs_idx_gq_agr = 1; // Separate index for AGR calculation

  for (j in 1:num_trees) {
    if (N_years_per_tree[j] > 1) {
      for (t_idx_local in 2:N_years_per_tree[j]) {
        real dt = years_per_tree[j, t_idx_local] - years_per_tree[j, t_idx_local - 1]; // Corrected years_per_tree here
        all_annual_agr_true[current_idx_agr] = (predicted_dbh_true_all[current_obs_idx_gq_agr + t_idx_local - 1] -
                                                predicted_dbh_true_all[current_obs_idx_gq_agr + t_idx_local - 2]) / dt;
        current_idx_agr += 1;
      }
    }
    current_obs_idx_gq_agr += N_years_per_tree[j];
  }

  // We can also calculate a "representative" mean AGR for the population
  real overall_mean_annual_agr = mean(AGR_drift_per_tree); // Mean of individual drifts
}

