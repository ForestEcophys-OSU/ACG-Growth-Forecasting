#' Estimate Tree DBH from Window Dendrometer Bands
#'
#' This function converts raw linear gap measurements from window-style manual dendrometer
#' bands into true cylindrical tree diameter (DBH) values. It utilizes numerical 
#' optimization to account for geometric surface curvature changes over time, following
#' the geometric principles used in the CTFS-ForestGEO network tracking protocols.
#'
#' @param dbh1 A numeric vector of initial or previous year tree diameters (in mm).
#' @param m1 A numeric vector of the previous year's dendrometer gap readings (in mm).
#' @param m2 A numeric vector of the current year's dendrometer gap readings (in mm).
#'
#' @return A numeric vector of estimated current tree diameters (in mm) of the same length 
#'   as the longest input vector. Returns \code{NA} for records with missing or invalid inputs.
#'
#' @export
#' @importFrom stats optimize
#'
#' @examples
#' # Single observation check
#' estimate_dendro_dbh(dbh1 = 72.0, m1 = 0.0, m2 = 3.8)
#'
#' # Vectorized processing check (e.g., within a mutate pipeline)
#' prev_dbhs <- c(72.0, 75.3, 79.1)
#' gaps_m1   <- c(0.0,  3.8,  6.2)
#' gaps_m2   <- c(3.8,  6.2,  9.9)
#' estimate_dendro_dbh(prev_dbhs, gaps_m1, gaps_m2)
estimate_dendro_dbh <- function(dbh1, m1, m2) {
  
  # =========================================================================
  # 1. INTERNAL GEOMETRIC AUXILIARY FUNCTIONS
  # =========================================================================
  
  # Objective/Loss function: evaluates how well a guessed diameter balances 
  # the theoretical arc-to-chord expansion equations against physical observations.
  objectiveFuncDendro <- function(diameter2, diameter1, gap1, gap2) {
    # Check for geometric impossibilities (a caliper gap cannot exceed the trunk diameter)
    if (gap1 > diameter1 || gap2 > diameter2) return(20)
    
    # Core mathematical identity: Conservation of band length over a cylinder surface.
    # Evaluates the difference between true cylindrical expansion and flat chord changes.
    delta <- abs(diameter1 - diameter2 + (1 / pi) * diameter2 * asin(gap2 / diameter2) - (1 / pi) * diameter1 * asin(gap1 / diameter1))
    return(delta)
  }
  
  # Single-record optimization solver: searches for the exact diameter that 
  # drives the objective function error down to its absolute minimum.
  findOneDendroDBH <- function(single_dbh1, single_m1, single_m2) {
    # Guard rail: return NA immediately if any variable required for the step is missing
    if (is.na(single_dbh1) || is.na(single_m1) || is.na(single_m2) || single_dbh1 <= 0) {
      return(NA_real_)
    }
    
    # Establish dynamic, bounded search boundaries to protect Brent's optimization method
    if (single_m2 > 0) {
      upper <- single_dbh1 + single_m2
    } else {
      upper <- single_dbh1 + 1
    }
    
    if (single_m2 < single_m1) {
      lower <- 0
    } else {
      lower <- single_dbh1
    }
    
    # Execute R's 1D numerical solver
    result <- stats::optimize(f = objectiveFuncDendro, 
                              interval = c(lower, upper), 
                              diameter1 = single_dbh1, 
                              gap1 = single_m1, 
                              gap2 = single_m2)
    return(result$minimum)
  }
  
  # =========================================================================
  # 2. VECTOR ALIGNMENT AND INPUT VALIDATION
  # =========================================================================
  
  # Find the longest vector to establish expected dataset length
  records <- max(length(dbh1), length(m1), length(m2))
  
  # HIGH-PERFORMANCE SHORTCUT: If inputs are single values, run instantly
  if (records == 1) {
    # Strip any indexing clutter and process immediately
    return(findOneDendroDBH(as.numeric(dbh1), as.numeric(m1), as.numeric(m2)))
  }
  
  # Standard R recycling rules: if an entry is length 1, clone it across the array
  if (length(dbh1) == 1) dbh1 <- rep(dbh1, records)
  if (length(m1) == 1)   m1   <- rep(m1,   records)
  if (length(m2) == 1)   m2   <- rep(m2,   records)
  
  # =========================================================================
  # 3. ROW-BY-ROW OPTIMIZATION EXECUTION (For true vectors)
  # =========================================================================
  
  # Pre-allocate an empty numeric output array for speed
  dbh2 <- numeric(records)
  
  # Step through the vector records and process each geometric puzzle independently
  for (i in 1:records) {
    dbh2[i] <- findOneDendroDBH(dbh1[i], m1[i], m2[i])
  }
  
  return(dbh2)
}
