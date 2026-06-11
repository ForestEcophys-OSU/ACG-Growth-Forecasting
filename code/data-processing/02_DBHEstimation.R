#'
#' # Dry Forest Growth Project
#' 
#' ## Script 1: Dendrometer band error calculation
#'
#' ### Author: German Vargas G.
#' 
#' ### Date: 2026-05-18
#' 
# Load packages ------
library(tidyverse)
library(ggpubr)

# Load functions ----
source("code/functions/calculate_DBH.R")# Function to calculate DBH from dendrometer data
source("code/functions/DBH_optimize.R")# Function to calculate DBH from dendrometer data

# Load datasets ------
# dendrometer database with linear measurements of 2024
tree.dbh <- read_csv("data/20260626_dendrometer_data_flagged.csv")

# Linear measurement and caliper comparison data collected in 2024
lm.correction.file <- read_csv(file = "data/LM_caliper_comparison_2024.csv") %>%# Load the linear measurement and caliper comparison data
  mutate(obsYear = rep(2024, nrow(.))) %>%# Add obsYear column with value 2024
  rename(linear_meas_mm = linear_measurement_mm, dendroRead_mm = caliper_mm,Sp_code.ORIG = Sp_code, DBH_tape_cm = diameter_2024_cm) %>%# Rename columns for consistency
  select(plot,subplot,tree,stem,dendrometer, Sp_code.ORIG, obsYear, DBH_tape_cm, dendroRead_mm, linear_meas_mm)# Select only the relevant columns for the correction file

# Add linear measurements of 2024 to the tree.dbh dataset by matching on plot, subplot, tree, stem, dendrometer and obsYear
tree.dbh <- tree.dbh %>%
  # 1. Bring in the linear_meas_mm column from the correction file
  left_join(lm.correction.file %>% select(plot, subplot, tree, stem, dendrometer, obsYear, lm_corr = linear_meas_mm),
            by = c("plot", "subplot", "tree", "stem", "dendrometer", "obsYear")) %>%
  # 2. If a match was found (lm_corr is not NA), overwrite the old value.
  #    Otherwise, keep the original linear_meas_mm value.
  mutate(linear_meas_mm = coalesce(lm_corr, linear_meas_mm)) %>%
  # 3. Clean up by removing the temporary join column
  select(-lm_corr)

rm(lm.correction.file) # Clean up the environment by removing the correction file)

# Calculate dbh for bias correction ------
# First, ensure your data is sorted chronologically per tree/stem
tree.dbh.final <- tree.dbh %>%
  # 1. Isolate target plot for assessing calculations
  #filter(plot == 1) %>%
  
  # 2. Explicitly arrange observations chronologically within every tree
  arrange(plot, tree, dendrometer, obsYear) %>%
  
  # 4. Construct the standard tidyverse indexing ID
  mutate(ID = paste(plot, tree, dendrometer, sep = "-")) %>%
  
  # 5. Group by tree ID and map our state engine safely over each split chunk
  group_by(ID) %>%
  group_modify(~ process_tree_dendro(.x, .y)) %>%
  
  # 6. Calculate the arc length transformation for all records where we have a valid dbh_dendro_mm and dendroRead_mm
  mutate(arclen_mm = acos(x = (1-((dendroRead_mm^2)/(2*((dbh_dendro_mm/2))^2))))*((dbh_dendro_mm/2))) %>%
  
  # 7. Ungroup and arrange the final dataset for sanity checks and visualizations
  ungroup() %>% arrange(plot, subplot, tree, dendrometer, obsYear)

# Assess the bias of dendrometer measurements -----
# number of linear measurements available for comparison
tree.dbh.final %>% filter(!is.na(linear_meas_mm) & !is.na(arclen_mm)) %>% nrow() # There are 273 records with linear measurements available for comparison.

summary(lm(tree.dbh.final$dendroRead_mm ~ tree.dbh.final$linear_meas_mm)) # The raw dendrometer readings show a strong correlation with the linear measurements, but the R^2 is slightly lower and the beta coefficient is less than 1, indicating that the dendrometer readings tend to underestimate DBH compared to the linear measurements. This suggests that there may be a bias in the raw dendrometer readings that causes them to consistently report smaller values than the linear measurements.
summary(lm(tree.dbh.final$arclen_mm~tree.dbh.final$linear_meas_mm))# The arc length transformation appears to have a slightly better fit (higher R^2 and beta closer to 1) compared to the raw dendrometer reading when regressed against the linear measurement. This suggests that the arc length transformation may be more effective at correcting for the bias in the dendrometer readings, providing a more accurate estimate of DBH that aligns closely with the linear measurements.

png(filename = "output/Fig_S_DendroReading_estimation.png", width = 10, height = 5, units = "in", res = 300)
ggarrange(
  tree.dbh.final %>% 
    ggplot(aes(y = dendroRead_mm, x = linear_meas_mm)) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_point(size=1.75,alpha=0.75) +
    stat_smooth(method = "lm",formula = "y~x", color = "blue") +
    annotate(geom="text", x=75, y=225, label=bquote(R^2==0.96*";"~beta==0.94),size=6)+
    xlim(0,260) + ylim(0,260) +
    labs(x = "Dendrometer reading (mm)", y = "Linear Measurement (mm)") +
    theme_bw() + # Start with a black and white theme for a cleaner look
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size=20), # Center and bold title
          panel.grid.major = element_blank(), # Remove major gridlines
          panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
          axis.line = element_line(color = "black"), # Add axis lines
          panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
          axis.text = element_text(size = 18), # Adjust axis text size
          axis.title = element_text(size = 20),
          strip.text = element_text(size=18),strip.background = element_blank()),
  tree.dbh.final %>% 
    ggplot(aes(y = arclen_mm, x = linear_meas_mm)) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_point(size=1.75,alpha=0.75) +
    stat_smooth(method = "lm",formula = "y~x", color = "blue") +
    annotate(geom="text", x=75, y=225, label=bquote(R^2==0.98*";"~beta==1.02),size=6)+
    xlim(0,260) + ylim(0,260) +
    labs(x = "Arc length (mm)", y = "Linear Measurement (mm)") +
    theme_bw() + # Start with a black and white theme for a cleaner look
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size=20), # Center and bold title
          panel.grid.major = element_blank(), # Remove major gridlines
          panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
          axis.line = element_line(color = "black"), # Add axis lines
          panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
          axis.text = element_text(size = 18), # Adjust axis text size
          axis.title = element_text(size = 20),
          strip.text = element_text(size=18),strip.background = element_blank()),
  labels = c("a)","b)"),ncol = 2, nrow = 1,font.label = list(size = 20))
dev.off()

# final estimation of DBH by year for all trees ------
tree.dbh.final <- tree.dbh.final %>%
  # 1. Rename variables to reflect the new meaning after transformation 
  # arclen_mm is now the corrected dendroRead_mm
  rename(dendroRead_mm.OLD = dendroRead_mm, dendroRead_mm = arclen_mm) %>%
  
  # 2. Explicitly arrange observations chronologically within every tree
  arrange(plot, tree, dendrometer, obsYear) %>%
  
  # 3. If the original dendroRead_mm is NA but we have a valid linear measurement, use the linear measurement 
  # as the new dendroRead_mm. Otherwise, keep the arclength based dendroRead_mm value.
  mutate(dendroRead_mm = if_else(condition = is.na(dendroRead_mm.OLD) & !is.na(linear_meas_mm), true = linear_meas_mm, false = dendroRead_mm)) %>% 
  
  # 4. Clean up by removing the previously estimated dbh_dendro_mm column, which is no longer needed after the arc length transformation
  select(-dbh_dendro_mm) %>%
  
  # 3. Group by tree ID and map our state engine safely over each split chunk
  group_by(ID) %>%
  group_modify(~ process_tree_dendro(.x, .y)) %>%
  mutate(dbh_dendro_cm = dbh_dendro_mm/10) %>%
  ungroup()

# visualize the final results for a sanity check for a randomly selected individual
random_tree <- tree.dbh.final %>% filter(ID == sample(unique(ID), 1))
ggplot(random_tree, aes(x = obsYear,y=dbh_dendro_mm/10)) +
  geom_line() +
  geom_point() +
  labs(x = "Year", y = "DBH (cm)",title = paste("Individual:",unique(random_tree$ID))) +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size=20), # Center and bold title
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 18), # Adjust axis text size
        axis.title = element_text(size = 20),
        strip.text = element_text(size=18),strip.background = element_blank())
rm(random_tree)

# save final results of DBH estimation by year -----
colnames(tree.dbh.final) # Check the column names of the final dataset to ensure they are correct and consistent with expectations.
tree.dbh.final %>%
  select(plot, subplot, tree, stem, dendrometer, ID, Sp_code, obsYear, dbh_dendro_cm,liana,phenology,canopy,mortality,QC,Notes) %>%
  write_csv("data/20260626_dbh_data_corrected.csv")

# datum <- read.csv("data/20260626_dbh_data_corrected.csv") # Read the saved CSV file to verify that it was saved correctly and the contents are as expecte
# datum %>% filter(obsYear >= 2021) %>%
#   select(plot, subplot, tree, stem, dendrometer, Sp_code, obsYear, dbh_dendro_cm,liana,phenology,canopy,mortality,QC,Notes) %>%
#   write_csv("public-data-releases/20260626_dbh_data_corrected.csv")