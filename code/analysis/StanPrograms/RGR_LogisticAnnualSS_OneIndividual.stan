// State-Space Tree Growth Model
//
// This model estimates latent daily/annual tree diameter (DBH) and relative 
// growth rate (RGR) from sequential dendrometer band increments. It explicitly 
// separates process variance (ecological growth dynamics) from observation error 
// (measurement noise).

data {
  int<lower=2> N;                // Total number of sequential observation steps
  real<lower=0> initial_tape_cm; // Initial diameter from manual tape (anchors the state)
  vector[N] y_delta;             // Observed dendrometer band increments (mm) at time t
}

parameters {
  // --- Global Population Parameters ---
  real mu_rgr;                   // Intrinsic base relative growth rate
  real beta_size;                // Density-dependent size effect on growth rate
  real<lower=0> sigma_proc;      // Standard deviation of process (ecological) noise
  real<lower=0> sigma_obs;       // Standard deviation of observation (measurement) error
  
  // --- Latent State Time-Series Variables ---
  vector[N] log_dbh;             // Unobserved true log-diameter (cm) at time t
  vector[N] rgr;                 // Realized relative growth rate at time t
}

transformed parameters {
  vector[N] predicted_delta_mm;  // Expected change in band reading (mm) between t-1 and t
  
  // Geometric correction: converts latent DBH (cm) differences to circumference increments (mm)
  // Circumference = DBH * pi. Factor of 10 converts cm to mm.
  for (t in 2:N) {
    predicted_delta_mm[t] = (exp(log_dbh[t]) - exp(log_dbh[t - 1])) * 3.1415926535 * 10;
  }
}

model {
  // --- 1. Initialization Block ---
  // Constrain the first latent state to be closely distributed around the manual tape measurement
  log_dbh[1] ~ normal(log(initial_tape_cm), 0.1); 
  
  // --- 2. Process Model (Latent State Evolution) ---
  for (t in 2:N) {
    // True growth rate is conditioned on the previous time step's size
    rgr[t] ~ normal(mu_rgr + beta_size * log_dbh[t - 1], sigma_proc);
    
    // Deterministic state update: new size = old size + growth rate
    // Modeled using a tight normal variance (1e-5) to satisfy Stan's parameter constraints
    log_dbh[t] ~ normal(log_dbh[t - 1] + rgr[t], 1e-5); 
  }
  
  // --- 3. Observation Model (Likelihood) ---
  // Evaluates observed dendrometer band values against the geometrically corrected latent growth
  y_delta[2:N] ~ normal(predicted_delta_mm[2:N], sigma_obs);
  
  // --- 4. Prior Distributions ---
  // Weakly informative structural priors for structural growth parameters
  mu_rgr ~ normal(0, 10);        
  beta_size ~ normal(0, 10);
  
  // Weakly informative scale-appropriate regularizing priors for variance parameters
  sigma_proc ~ exponential(1);   
  sigma_obs ~ exponential(1);    
}