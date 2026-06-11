# Load libraries
library(cmdstanr)
library(tidyverse)

# 1. Filter data for a single individual 
# Ensure it is ordered sequentially by time step (year/day)
one_tree_data <- tree_data %>%
  filter(tree_id == "Tree_01") %>%  # Change to your identifier
  arrange(time_step)

# 2. Format the data into a list for Stan
stan_data <- list(
  N               = nrow(one_tree_data),
  initial_tape_cm = one_tree_data$DBH_tape_cm[1], # First manual measurement
  y_delta         = one_tree_data$delta_dendro_mm # Vector of increments
)

# 3. Compile the Stan model
# This translates the Stan file to C++ and compiles it
model <- cmdstan_model("")

# 4. Sample from the posterior
fit <- model$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 200 # Print updates every 200 iterations
)