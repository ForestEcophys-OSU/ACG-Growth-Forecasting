data {
  int<lower=1> N; // Total number of diameter measurements across all trees
  int<lower=1> N_trees; // Total number of individual trees
  int<lower=1> tree_N_obs[N_trees]; // Number of observations for each tree
  int<lower=1, upper=N> tree_start_idxs[N_trees]; // Starting index for each tree's data in the flat vectors
  int<lower=1, upper=N> tree_end_idxs[N_trees]; // Ending index for each tree's data
  vector[N] tree_years; // Concatenated vector of years for all tree observations
  vector[N] tree_dbhs; // Concatenated vector of DBH for all tree observations
  
  int<lower=1> N_year_grid; // Number of points in the prediction year grid
  vector[N_year_grid] year_grid; // The year grid for predictions
}

transformed data {
  // Pre-calculate natural logarithm of all DBH measurements for efficiency
  vector[N] log_tree_dbhs = log(tree_dbhs); 
}

parameters {
  //real<lower=0> d0[N_trees]; // Diameter at year of first observation (cm)
  //real beta[N_trees];        // Linear growth rate (cm / year)
  real r1;         // RGR parameter 1 (intercept for log-initial DBH)
  real r2;         // RGR parameter 2 (slope for log-initial DBH)
  real<lower=0> sigma;       // Measurement variability (cm) - additive error on original DBH scale
}

model {
  // Priors for r1 and r2
  // Assuming variance of 10^4 (as per "N(0, 10-4)" and ignoring the 10^-4 notation ambiguity,
  // interpreting it as a large variance for non-informative prior), so std dev of 100.
  r1 ~ normal(0, 100); 
  r2 ~ normal(0, 100); 
  
  //d0 ~ normal(0, 150 / 2.57);     // 99% prior mass between 0 and 150 cm
  //beta ~ normal(0, 2 / 2.32);     // 99% prior mass between +/- 2 cm / year
  
  // Prior for sigma: 99% prior mass between 0 and 0.25 cm
  // This prior implies sigma is an additive error on the original DBH scale.
  sigma ~ normal(0, 0.25 / 2.57); 
  
  for (t in 1:N_trees) {
    int start_idx = tree_start_idxs[t];
    int end_idx = tree_end_idxs[t];
    
    // Extract observations for the current tree
    vector[tree_N_obs[t]] years = tree_years[start_idx:end_idx];
    vector[tree_N_obs[t]] log_dbhs = log_tree_dbhs[start_idx:end_idx];
    vector[tree_N_obs[t]] current_tree_dbhs = tree_dbhs[start_idx:end_idx]; // Original DBH for likelihood
    //vector[tree_N_obs[t]] dbhs = tree_dbhs[start_idx:end_idx];
    
    // Initial diameter for tree t at its first observation
    real log_D1i = log_dbhs[1];
    
    // Calculate the Relative Growth Rate (Ri) for this specific tree
    // Ri = r1 + r2 * ln(D_1i) as in Iida et al. 2014
    real Ri = r1 + r2 * log_D1i;
    
    // Calculate the predicted logarithm of DBH for each observation of this tree
    // ln(D_2i) = ln(D_1i) + Ri * (t_2i - t_1i)
    vector[tree_N_obs[t]] log_dbh_predicted;
    
    for (i in 1:tree_N_obs[t]) {
      log_dbh_predicted[i] = log_D1i + Ri * (years[i] - years[1]);
    }
    
    // Likelihood:
    // The observed DBHs (on the original scale) are normally distributed
    // around the exponentiated predicted log-DBHs (also on original scale),
    // with a standard deviation 'sigma' (measurement error in cm).
    current_tree_dbhs ~ normal(exp(log_dbh_predicted), sigma);
    #dbhs ~ normal(d0[t] + beta[t] * (years - years[1]), sigma);
  }
}

generated quantities {
  // This block is for generating predictions or other quantities of interest
  // after the model has been fit.
  // For example, to predict DBH on the 'year_grid':
  // You would typically define a starting DBH (e.g., an average initial DBH)
  // for which you want to make predictions.
  
  // Declare the array for predicted DBHs on the grid
  real dbh_grid_pred[N_trees, N_year_grid];
  
  // Loop through each tree to make predictions
  for (t in 1:N_trees) {
    // Get the initial year and log(DBH) for the current tree
    real y0 = tree_years[tree_start_idxs[t]];
    real log_D1i_current_tree = log_tree_dbhs[tree_start_idxs[t]]; // This is log(D1i) for tree t
    
    // Calculate the RGR (Ri) for the current tree based on its initial log(DBH)
    real Ri_current_tree = r1 + r2 * log_D1i_current_tree;
    
    for (n in 1:N_year_grid) {
      // Calculate the mean of log(DBH) at year_grid[n] for the current tree
      // ln(D_grid) = ln(D_1i) + Ri * (t_grid - t_1i)
      real log_mu_pred = log_D1i_current_tree + Ri_current_tree * (year_grid[n] - y0);
      
      // Convert the mean of log(DBH) back to the original DBH scale
      real mu_pred = exp(log_mu_pred);
      
      // Sample a predicted DBH value from the normal distribution
      // This incorporates the 'sigma' measurement error on the original DBH scale
      dbh_grid_pred[t, n] = normal_rng(mu_pred, sigma);
    }
  }
}
