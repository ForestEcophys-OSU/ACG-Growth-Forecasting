// Stan code for Bayesian State-Space Tree Growth Model
// This model infers latent annual increments and true diameters
// based solely on periodic diameter measurements.

data {
  int<lower = 0> N_trees;               // Total number of unique trees
  int<lower = 0> N_years;               // Total number of years in the study period
  int<lower = 0> N_dbh;                 // Total number of dbh observations from dendrometer bands
  
  // Indices and values for diameter measurements
  int<lower = 0> tree_dbh[N_dbh]; // Index of dbh measurements denoting the tree number 
  int<lower = 0> year_dbh[N_dbh]; // Index of dbh measurements denoting the year number
  real<lower = 0> dbh[N_dbh]; // dbh measurement vector
}

parameters {
  // Variance parameters
  real<lower=0, upper =10> sigma_dbh; // Measurement error for diameter observations
  
  // Tree-level random effect parameters for the mean of log-increments
  real<lower=-10, upper =10> mutree;     // Mean of the tree-level random effects
  real alpha_TREE[N_trees];           // Tree-specific random effects for the mean of log-increments
  real<lower=1e-6> sigma_TREE;           // Standard deviation of the tree-level random effects
  
  // Latent true initial diameter for each tree (at year 0 or conceptual start)
  real<lower=0> dbh_init[N_trees];
}

transformed parameters {
  // True diameter 'x' is a transformed parameter, derived from xinit and increments
  // This matrix will hold the estimated true diameter for each tree in each year.
  matrix<lower=1e-6> DBH[N_trees, N_years];

  for (i in 1:N_trees) {
    // Diameter at year 1 is initial diameter plus the first year's increment
    DBH[i,1] = DBHinit[i] + inc[i, 1];

    // Subsequent diameters are cumulative sums of increments
    for (t in 2:N_years) {
      DBH[i, t] = DBH[i, t-1] + inc[i, t];
    }
  }
}

transformed parameters{
   real<lower=1e-6> x[Nrow, Ncol]; //true diameter x is a transformed parameter

  for (i in 1:Nrow){

    x[i,1] = xinit[i] + inc[i, 1];

    for (t in 2:Ncol) {
      x[i, t] = x[i, t-1] + inc[i, t] ;
    }
}

model {
  // Priors on tree-level random effect (for the mean of log-increments)
  mutree ~ normal(0, 5);
  sigma_TREE ~ cauchy(0, 5);
  alpha_TREE ~ normal(mutree, sigma_TREE);
  
  // Variance priors
  // sigma_dbh: Measurement error for diameter observations
  sigma_dbh ~ normal(1, 0.01); // Centered at 1cm, very tight prior. Adjust if your measurement error is smaller.
                              // (e.g., normal(0.1, 0.01) if you expect ~0.1cm error)

  // sigma_growth: Process error for the true annual increments
  sigma_growth ~ normal(0, 5); // Wide prior, allowing for substantial variability in increments
  
  // xinit initial prior - weakly informative
  // This prior allows Stan to estimate the initial true diameter for each tree.
  DBHinit ~ uniform(0, 200); // Uniform prior between 0 and 200 cm for initial diameter
  
  // Tree diameter increment process model
  // This defines how the true annual increments (`inc`) evolve.
  // Lognormal ensures increments are positive, `alpha_TREE` provides tree-specific mean.
  for (i in 1:N_trees) {
    inc[i,1:N_years] ~ lognormal(alpha_TREE[i], sigma_growth);
  }
  
  // Diameter data model
  // This links the observed diameter measurements (`z`) to the derived true diameters (`x`).
  for (d in 1:N_dbh) {
    dbh[d] ~ normal(DBHinit[tree_dbh[d], year_dbh[d]], sigma_dbh);
  }
}

   for(d in 1:Ndia){
     z[d] ~ normal(x[treedia[d], yeardia[d]], sigma_dbh);
     //diameter data z corresponds to the true diameter matrix through the tree and year diameter indices
   }





}
model {
  
  //priors on tree-level random effect
  mutree ~ normal(0, 5);
  sigma_TREE ~ cauchy(0, 5);
  alpha_TREE ~ normal(mutree, sigma_TREE);
  
  
  //variance priors
  sigma_dbh ~ normal(1, 0.01); //normal(1,0.01) works well on base model
  sigma_inc ~ normal(0.035, 0.01); //based on SD of remeasured increments
  sigma_add ~ normal(0, 5); // wide (ish) prior given that max diameter inc is ~3.5
  
  //x initial prior-weakly informative --max tree size in data is 75
  xinit  ~ uniform(0, 75);
 
  
 // tree ring diameter increment process model
 for(i in 1:Nrow){//loop over the number of trees
           //} changed dec 20
       inc[i,1:Ncol] ~ lognormal( alpha_TREE[i] , sigma_add);
    
  }
  
  //tree ring increment data model
  for(j in 1:Ninc){ //loop over the number of increment measurements
     y[j] ~ normal(inc[treeinc[j], yearinc[j]], sigma_inc)T[0,];
     //increment data y corresponds to the true increment matrix through the tree and year increment indices
   }
 
 //diameter data model --moved from previous loop
   //  z[i,t] ~ normal(x[i,t], sigma_dbh);
   for(d in 1:Ndia){
     z[d] ~ normal(x[treedia[d], yeardia[d]], sigma_dbh);
     //diameter data z corresponds to the true diameter matrix through the tree and year diameter indices
   }
}