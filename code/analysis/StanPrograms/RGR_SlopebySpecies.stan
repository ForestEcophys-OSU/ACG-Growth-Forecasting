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

  // NEW DATA: Species information - REQUIRED FOR SPECIES-LEVEL HIERARCHY
  int<lower=1> N_species; // Total number of species
  int<lower=1, upper=N_species> tree_species_idx[N_trees]; // Index of the species for each tree
}

transformed data {
  // Pre-calculate natural logarithm of all DBH measurements for efficiency
  vector[N] log_tree_dbhs = log(tree_dbhs); 
}

parameters {
  // Raw parameters for non-centered parameterization for *species-level* r1 and r2
  // These replace the global 'r1' and 'r2' from the equal slope model
  vector[N_species] r1_species_raw; 
  vector[N_species] r2_species_raw; 
  
  real<lower=0> sigma; // Measurement variability (cm)
  
  // Hyperparameters for the *community-level* mean and variability of species r1 and r2
  real mu_r1_community;
  real<lower=0> sigma_r1_species; // Standard deviation of r1 across species
  real mu_r2_community;
  real<lower=0> sigma_r2_species; // Standard deviation of r2 across species
}

transformed parameters {
  // Species-specific RGR parameters (r1j, r2j in Iida et al. 2014)
  // These are derived from the community-level means and species-level raw deviations
  real r1_species[N_species]; 
  real r2_species[N_species]; 
  
  for (j in 1:N_species) {
    r1_species[j] = mu_r1_community + sigma_r1_species * r1_species_raw[j]; // Non-centered reparameterization
    r2_species[j] = mu_r2_community + sigma_r2_species * r2_species_raw[j]; // Non-centered reparameterization
  }
}

model {
  // Priors for community-level hyperparameters (mu_r1_community, mu_r2_community)
  // Iida et al. (2014) use N(0, 10^4) for these. Using standard deviation of 100.
  mu_r1_community ~ normal(0, 100); 
  mu_r2_community ~ normal(0, 100);
  
  // Priors for species-level standard deviations (sigma_r1_species, sigma_r2_species)
  // Iida et al. (2014) use Gamma(10^-2, 10^-2) for variance. For standard deviation,
  // a Half-Normal or Half-Cauchy is often preferred in Stan for better sampling.
  // Sticking with normal(0, 10) for now, constrained positive by 'lower=0'.
  sigma_r1_species ~ normal(0, 10); 
  sigma_r2_species ~ normal(0, 10);
  
  // Priors for the raw species parameters (standard normal, as per non-centered reparameterization)
  r1_species_raw ~ normal(0, 1);
  r2_species_raw ~ normal(0, 1);
  
  // Prior for sigma (measurement error)
  sigma ~ normal(0, 0.25 / 2.57); 
  
  for (t in 1:N_trees) {
    int start_idx = tree_start_idxs[t];
    int end_idx = tree_end_idxs[t];
    
    // Get the species index for the current tree
    int current_species_idx = tree_species_idx[t];

    vector[tree_N_obs[t]] years = tree_years[start_idx:end_idx];
    vector[tree_N_obs[t]] log_dbhs = log_tree_dbhs[start_idx:end_idx];
    vector[tree_N_obs[t]] current_tree_dbhs = tree_dbhs[start_idx:end_idx]; 
    
    real log_D1i = log_dbhs[1]; // Initial diameter for tree t at its first observation
    
    // Calculate the Relative Growth Rate (Ri) for this specific tree based on its species' parameters
    // Ri = r1j + r2j * ln(D1i) as in Iida et al. (2014)
    real Ri = r1_species[current_species_idx] + r2_species[current_species_idx] * log_D1i;
    
    // Calculate the predicted logarithm of DBH for each observation of this tree
    // ln(D2i) = ln(D1i) + Ri * (t2i - t1i)
    vector[tree_N_obs[t]] log_dbh_predicted;
    for (i in 1:tree_N_obs[t]) {
      log_dbh_predicted[i] = log_D1i + Ri * (years[i] - years[1]);
    }
    
    // Likelihood: Observed DBHs (on original scale) are normally distributed
    // around the exponentiated predicted log-DBHs, with standard deviation 'sigma'.
    current_tree_dbhs ~ normal(exp(log_dbh_predicted), sigma);
  }
}

generated quantities {
  real dbh_grid_pred[N_trees, N_year_grid];
  
  for (t in 1:N_trees) {
    int current_species_idx = tree_species_idx[t]; // Get species index for this tree
    real y0 = tree_years[tree_start_idxs[t]];
    real log_D1i_current_tree = log_tree_dbhs[tree_start_idxs[t]]; 
    
    // Calculate the RGR (Ri) for the current tree on the prediction grid, based on its species' parameters
    real Ri_current_tree = r1_species[current_species_idx] + r2_species[current_species_idx] * log_D1i_current_tree;
    
    for (n in 1:N_year_grid) {
      real current_year_grid = year_grid[n];
      
      real log_mu_pred = log_D1i_current_tree + Ri_current_tree * (current_year_grid - y0);
      real mu_pred = exp(log_mu_pred);
      
      dbh_grid_pred[t, n] = normal_rng(mu_pred, sigma);
    }
  }
}
