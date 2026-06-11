#' Calculate Rainfall Seasonality and Drought Metrics
#'
#' This function processes a matrix of monthly rainfall data to align it by the 
#' site's specific biological/hydrological year rather than the calendar year. It then 
#' computes several advanced climatic indices including Feng's entropy/seasonality metrics, 
#' and Aragao's Maximum Climatological Water Deficit (MCWD).
#'
#' @param data A data frame or matrix where the first column is the calendar year, 
#'   and columns 2 to 13 are the 12 calendar months (Jan to Dec).
#' @param Site A character string specifying the name of the study site or station.
#'
#' @return A data frame containing annual metrics for each calculated hydrological year:
#' \item{site}{Name of the site}
#' \item{year}{The starting calendar year of the hydrological period}
#' \item{hydro_start_mos}{The calendar month that marks the start of this hydrological year}
#' \item{annual_rainfall}{Total rainfall accumulated in the 12-month hydrological year (MAR)}
#' \item{entropy}{Relative entropy (D), measuring the concentration of precipitation (Feng et al. 2013)}
#' \item{seasonality}{Seasonality index (S), scaling entropy by overall relative abundance}
#' \item{spread}{The standard deviation/duration of the rainy season around the centroid}
#' \item{centroid}{The average timing/center of mass of the rainy season (month 1 to 12)}
#' \item{interwet_months}{The number of distinct wet periods separated by dry spells (modality)}
#' \item{mcwd}{Maximum Climatological Water Deficit (mm), tracking cumulative water stress (Aragao et al. 2007)}
#' \item{drymonths}{The count of months during the hydrological cycle tracking a water deficit}
#' 
#' @importFrom zoo rollapply
#' @export
seasonality <- function(data = NULL, Site) {
  
  # =========================================================================
  # 1. HYDROLOGICAL YEAR ALIGNMENT
  # =========================================================================
  
  # Identify the driest calendar month of each year to anchor the tracking frame.
  # The hydrological cycle is defined to start the month *after* the absolute minimum rainfall.
  hydroYearStart <- apply(data[, 2:13], 1, function(x) which(x == min(x))[1]) + 1
  
  years <- data[, 1] # Assuming the first column contains the year labels
  
  m <- as.matrix(data[, 2:13]) # Convert the monthly rainfall data to a matrix, excluding the year column
  #rownames(m) <- years # Set the row names to the corresponding years for easier tracking
  
  # Initialize an empty matrix to hold the shifted, true hydrological years.
  # Note: The final calendar year is dropped because its cycle bleeds into the missing next year.
  rainfall_mat <- matrix(nrow = (nrow(data) - 1), ncol = 12,
                         dimnames = list(years[1:(length(years) - 1)], c(1:12)))
  
  # Vectors to dynamically record tracking information for each cycle
  saved_start_years <- c()
  saved_start_months <- c()
  
  # Re-align monthly columns so that Month 1 is always the start of the site's wet build-up
  for(i in 1:(length(m[,1]) - 1)) {
    if(hydroYearStart[i] != 13) {
      a <- 13 - hydroYearStart[i]
      # Pull months from the current calendar year (from dry minimum to December)
      rainfall_mat[i, 1:a] <- m[i, hydroYearStart[i]:12]
      
      # Fill the remaining months using the beginning of the subsequent calendar year
      b <- a + 1
      c <- 12 - a
      rainfall_mat[i, b:12] <- m[i + 1, 1:c]
      
      saved_start_years[i] <- years[i]
      saved_start_months[i] <- hydroYearStart[i]+1 # The month immediately following the dry minimum
    } else {
      # If the minimum month was December, the next biological year maps perfectly to the next calendar year
      rainfall_mat[i, 1:12] <- m[i + 1, 1:12] 
      
      saved_start_years[i] <- years[i + 1]
      saved_start_months[i] <- 1 # January
    }
  }
  
  # =========================================================================
  # 2. FENG ET AL. (2013) SEASONALITY & ENTROPY METRICS
  # =========================================================================
  
  # Helper function to compute the directional standard deviation (spread) 
  # of rainfall distributed around the seasonal centroid.
  spreadFunc <- function(row) {
    return(sqrt((1 / AR[row]) * sum(((mos - cent[row])^2) * rainfall_mat[row, ])))
  }
  
  # AR: Mean Annual Rainfall (MAR) calculated across the customized hydrological footprint
  AR <- rowSums(rainfall_mat)
  
  # pm: Normalizes rainfall into monthly fractional percentages sum up to 1.0
  pm <- (rainfall_mat) / AR
  
  # D: Information Entropy / Relative Entropy. Measures how uniform or concentrated rain is.
  D <- rowSums(pm * log2(pm / (1 / 12) + 0.0001))
  
  # seas: Seasonality Index scaling structural entropy by relative precipitation abundance
  seas <- D * AR / max(AR)
  
  # cent: The "center of mass" or timing centroid of the rainy season (valued 1 to 12)
  cent <- (1 / AR) * apply(rainfall_mat, 1, function(x) sum(x * 1:12))
  mos <- 1:12
  
  # spr: Evaluates the duration/spread envelope of the wet season block
  spr <- sapply(1:nrow(rainfall_mat), spreadFunc)
  
  # =========================================================================
  # 3. MODALITY & WET/DRY TRANSITIONS
  # =========================================================================
  
  # Establish a logical grid identifying wet months (> 100mm/month standard threshold)
  wetMat <- rainfall_mat > 100 
  
  # nmodes: Calculates if a year is bi-modal or multi-modal by tracking transitions.
  nmodes <- apply(wetMat, 1, function(x) {
    sum(zoo::rollapply(as.numeric(rle(x)$values), width = 3, identical, y = c(1, 0, 1)))
  })
  
  # =========================================================================
  # 4. ARAGAO ET AL. (2007) MAXIMUM CLIMATOLOGICAL WATER DEFICIT (MCWD)
  # =========================================================================
  
  # An iterative tracking calculation evaluating cumulative hydrological stress.
  MCWDfunc <- function(rainRow) {
    waterDef <- NULL
    waterDef[1] <- 0  # Anchor initial boundary state
    
    for(j in 2:12) {
      if((waterDef[j - 1] - 100 + rainRow[j]) < 0) {
        waterDef[j] <- waterDef[j - 1] - 100 + rainRow[j]
      } else {
        waterDef[j] <- 0
      }
    }
    MCWD <- min(waterDef)          
    # Calculate the maximum CONSECUTIVE dry months using Run Length Encoding (rle)
    is_dry <- waterDef < 0
    if (any(is_dry)) {
      dry_runs <- rle(is_dry)
      # Extract lengths where the run value is TRUE
      max_consecutive_dry <- max(dry_runs$lengths[dry_runs$values == TRUE])
    } else {
      max_consecutive_dry <- 0
    }
    return(c(MCWD, max_consecutive_dry))
  }
  
  # Execute the custom water deficit tracking loop row-by-row across the aligned matrices
  yearMCWD <- t(apply(rainfall_mat, 1, MCWDfunc))
  
  # =========================================================================
  # 5. DATA COMPILING & EXPORT
  # =========================================================================
  
  # Translate numeric month coordinates to tidy abbreviations
  #month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  #saved_start_months_label <- month_names[saved_start_months]
  
  # Bundle all localized vectors cleanly into a standardized output layout
  Seasonality_metrics <- data.frame(
    site            = rep(Site, nrow(rainfall_mat)),
    year            = saved_start_years,             # Uses recovered vector instead of row.names
    hydro_start_mos = saved_start_months,      # Appends new month character tracking vector
    annual_rainfall = AR,
    #pm              = pm,                   # monthly fraction for reference
    entropy         = D,                             
    seasonality     = seas,
    spread          = spr,
    centroid        = cent,
    #interwet_months = nmodes,
    mcwd            = yearMCWD[, 1],
    #drymonths       = yearMCWD[, 2],
    stringsAsFactors = FALSE
  )
  
  return(Seasonality_metrics)
}

