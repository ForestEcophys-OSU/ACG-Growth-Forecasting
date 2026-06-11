library(dplyr)
library(purrr)

process_tree_dendro <- function(df, key) {
  # Pre-allocate output columns within this tree's context
  df$dbh_dendro_mm <- NA_real_
  #df$dbh_cm        <- NA_real_
  
  # Initialize state memory variables
  last_valid_idx <- NA
  last_valid_m   <- NA
  
  # Tracks what band we are currently working on
  current_band   <- df$dendroNumber[1]
  
  # Flag to detect when a band has shifted, but hasn't received its first valid value yet
  band_init_pending <- FALSE 
  
  # Loop over the rows of this individual tree
  for (i in 1:nrow(df)) {
    
    # Check if the row's band changed from our tracked current_band
    if (df$dendroNumber[i] != current_band) {
      current_band      <- df$dendroNumber[i]
      band_init_pending <- TRUE # A new band has started; its initialization is now pending
    }
    
    # -----------------------------------------------------------------------
    # INITIAL REGISTRATION: First row of this specific tree
    # -----------------------------------------------------------------------
    if (i == 1) {
      df$dbh_dendro_mm[i] <- df$DBH_tape_cm[i] * 10
      #df$dbh_cm[i]        <- df$DBH_tape_cm[i]
      df$dendroRead_mm[i] <- 0 
      
      last_valid_idx    <- i
      last_valid_m      <- 0
      band_init_pending <- FALSE # Row 1 acts as its own baseline initialization
    } 
    # -----------------------------------------------------------------------
    # SEQUENTIAL TRACKING: Same band or handling a delayed band initialization
    # -----------------------------------------------------------------------
    else {
      # Skip row if it's missing data or marked QC == -1
      if (is.na(df$dendroRead_mm[i]) || (!is.na(df$QC[i]) && df$QC[i] == -1)) {
        next 
      } 
      
      # ---------------------------------------------------------------------
      # BRANCH A: This is the FIRST valid measurement of a new band!
      # (Runs immediately even if previous years on this band number were NA)
      # ---------------------------------------------------------------------
      else if (band_init_pending) {
        base_dbh_mm <- df$dbh_dendro_mm[last_valid_idx]
        #base_dbh_cm <- df$dbh_cm[last_valid_idx]
        
        # New band initialization baseline always treats m1 as 0
        df$dbh_dendro_mm[i] <- estimate_dendro_dbh(
          dbh1 = base_dbh_mm, 
          m1   = 0, 
          m2   = df$dendroRead_mm[i]
        )
        
        #delta_inc_mm <- df$dendroRead_mm[i] - 0
        #df$dbh_cm[i] <- base_dbh_cm + (delta_inc_mm / (10 * pi))
        
        # Save state trackers and turn off the initialization pending flag
        last_valid_idx    <- i
        last_valid_m      <- df$dendroRead_mm[i]
        band_init_pending <- FALSE 
      } 
      
      # ---------------------------------------------------------------------
      # BRANCH B: Sequential tracking on an already initialized band
      # ---------------------------------------------------------------------
      else {
        current_m  <- df$dendroRead_mm[i]
        target_idx <- last_valid_idx
        target_m   <- last_valid_m
        
        # Handle Decreases via backward scanning
        if (current_m < target_m) {
          found_smaller <- FALSE
          
          for (j in (i-1):1) {
            if (df$dendroNumber[j] == df$dendroNumber[i] && 
                !is.na(df$dbh_dendro_mm[j]) && 
                !is.na(df$dendroRead_mm[j]) && 
                df$dendroRead_mm[j] < current_m) {
              
              target_idx    <- j
              target_m      <- df$dendroRead_mm[j]
              found_smaller <- TRUE
              break 
            }
          }
          
          if (!found_smaller) {
            # Find the true initialization row for this specific band
            band_start_idx <- which(df$dendroNumber == df$dendroNumber[i] & !is.na(df$dbh_dendro_mm))[1]
            target_idx     <- band_start_idx
            target_m       <- 0 
          }
        }
        
        # Run calculations
        df$dbh_dendro_mm[i] <- estimate_dendro_dbh(
          dbh1 = df$dbh_dendro_mm[target_idx], 
          m1   = target_m, 
          m2   = current_m
        )
        
        #delta_inc_mm <- current_m - target_m
        #df$dbh_cm[i] <- df$dbh_cm[target_idx] + (delta_inc_mm / (10 * pi))
        
        # Update running tracking metrics
        last_valid_idx <- i
        last_valid_m   <- current_m
      }
    }
  }
  
  return(df)
}