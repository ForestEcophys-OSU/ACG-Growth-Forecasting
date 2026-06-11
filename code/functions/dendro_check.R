#' dendro.check
#'
#' Interactive quality control check for dendrometer increment data. 
#' Allows the user to visually inspect and manually flag outliers based on dendrometer readings 
#'
#' @param datum A data frame containing the following columns:
#'   - `obsYear`: Year of observation
#'   - `dendroRead_mm`: Dendrometer readings (NA if not available)
#'   - `linear_meas_mm`: Linear measurements of dendrometer increment (NA if not available)
#'   - `dendrometer`: Unique identifier for the tree dendrometer band
#'   - `Sp_code`: Species code column
#'   - `dendroNumber`: Dendrometer band tracking number
#'   - `newDendrometer`: Dendrometer band version/reset indicator
#'   - `mortality`: Binary indicator (1 if tree died, 0 or NA otherwise)
#'   - `QC` (optional): Column indicating quality control status; if missing, will be created
#'
#' @return The same data frame as input, invisibly returned, with an updated `QC` column
#' @export
dendro.check <- function (datum) {
  
  # 1. Initialize Quality Check column (QC) if missing
  if (!("QC" %in% names(datum))) {
    datum$QC <- rep(0, nrow(datum))
  }
  
  # ==========================================================
  # AUTOMATED QC PRESETS (EFFICIENCY UPDATES)
  # ==========================================================
  
  # Preset A: Automatically flag known mortality events as Fail (-1)
  if ("mortality" %in% names(datum)) {
    mortality_rows <- which(datum$mortality == 1)
    if (length(mortality_rows) > 0) {
      datum$QC[mortality_rows] <- -1
    }
  }
  
  # Preset B: If there is only 1 unique year of data, auto-assign Pass (1) and skip plot
  unique_years <- unique(datum$obsYear[!is.na(datum$obsYear)])
  if (length(unique_years) <= 1) {
    # Set any rows that aren't already flagged as mortality to Pass (1)
    untouched_rows <- which(datum$QC == 0)
    if (length(untouched_rows) > 0) {
      datum$QC[untouched_rows] <- 1
    }
    return(invisible(datum))
  }
  
  # ==========================================================
  # STANDARDIZED DATA PREPARATION
  # ==========================================================
  
  # 2. Skip condition: Check if everything is NA or insufficient observations exist
  all_dendro_na <- all(is.na(datum$dendroRead_mm))
  all_linear_na <- all(is.na(datum$linear_meas_mm))
  
  if ((all_dendro_na && all_linear_na) || nrow(datum) <= 1) {
    datum$QC[datum$QC == 0] <- 1
    return(invisible(datum))
  }
  
  # Construct a dynamic title using species, band number, and version/reset status
  sp_code <- if("Sp_code" %in% names(datum)) datum$Sp_code[1] else "UnknownSp"
  d_num   <- if("dendrometer" %in% names(datum)) datum$dendrometer[1] else "NoNum"
  d_ver   <- if("dendroNumber" %in% names(datum)) datum$dendroNumber[1] else "V1"
  
  meta_title <- paste0("Species: ", sp_code, " | Band #: ", d_num, " | Version: ", d_ver)
  
  # UI Parameters
  status <- c("red", "black", "blue", "grey") # -1 = fail, 0 = unchecked, 1 = pass
  error_threshold <- 10 # mm threshold for visual error bars
  
  # Obtain last year in the census
  latest_year <- max(datum$obsYear, na.rm = TRUE)
  
  # 3. Build a coordinate map of valid interactive points across columns
  id_map <- data.frame(row_idx = integer(), value = numeric(), type = character())
  for (i in 1:nrow(datum)) {
    if (!is.na(datum$dendroRead_mm[i])) {
      id_map <- rbind(id_map, data.frame(row_idx = i, value = datum$dendroRead_mm[i], type = "dendro"))
    }
    if (!is.na(datum$linear_meas_mm[i])) {
      id_map <- rbind(id_map, data.frame(row_idx = i, value = datum$linear_meas_mm[i], type = "linear"))
    }
  }
  
  if (nrow(id_map) == 0) {
    datum$QC[datum$QC == 0] <- 1
    return(invisible(datum))
  }
  
  # Force strict axes limits to move the "FAIL ALL" button into the true upper left corner
  x_limits <- c(2008, latest_year)
  y_range  <- range(c(datum$dendroRead_mm, datum$linear_meas_mm), na.rm = TRUE)
  
  # Expand limits to ensure the visual error whiskers are fully contained inside the window bounds
  y_limits <- c(y_range[1] - error_threshold - 5, y_range[2] + error_threshold + 5)
  if(any(is.infinite(y_limits))) y_limits <- c(0, 100)
  
  # "FAIL ALL" Button Coordinates (True upper left corner of plot bounds)
  ulhc_x <- x_limits[1] + 0.5
  ulhc_y <- y_limits[2] - (0.05 * diff(y_limits))
  
  # Append button coordinates to our interactive registry pool
  click_x <- c(datum$obsYear[id_map$row_idx], ulhc_x)
  click_y <- c(id_map$value, ulhc_y)
  
  # ==========================================================
  # PLOT 1: SELECT OUTLIERS
  # ==========================================================
  
  # Print instructions directly to the R console
  cat("\n------------------------------------------------------------\n")
  cat("QC CONSOLE INSTRUCTIONS - STEP 1 (OUTLIER SELECTION):\n")
  cat("1. Look at the graphics window. Light gray bars show the +/- 10mm buffer.\n")
  cat("2. Click LEFT on individual points to flag them as OUTLIERS.\n")
  cat("   - Circles (●) represent standard band readings.\n")
  cat("   - Triangles (▲) represent manual linear measurements.\n")
  cat("3. Click the red 'FAIL ALL' text in the upper left corner to invalidate the entire series.\n")
  cat("4. When finished flagging points for this tree, press ESC or RIGHT-CLICK and select 'Stop'.\n")
  cat("------------------------------------------------------------\n\n")
  
  plot(NULL, NULL, 
       xlim = x_limits, ylim = y_limits,
       xlab = "Observation Year", ylab = "Increment Measurement (mm)",
       main = paste("CLICK ON OUTLIERS\n", meta_title))
  
  # Background Error Whiskers Strategy (Rendered under data points)
  if (!all_dendro_na) {
    valid_dendro <- !is.na(datum$dendroRead_mm)
    graphics::arrows(datum$obsYear[valid_dendro], datum$dendroRead_mm[valid_dendro] - error_threshold,
                     datum$obsYear[valid_dendro], datum$dendroRead_mm[valid_dendro] + error_threshold,
                     length = 0.05, angle = 90, code = 3, col = "gray80", lwd = 1.5)
    points(datum$obsYear, datum$dendroRead_mm, col = status[datum$QC + 2], pch = 20, cex = 2)
  }
  if (!all_linear_na) {
    valid_linear <- !is.na(datum$linear_meas_mm)
    graphics::arrows(datum$obsYear[valid_linear], datum$linear_meas_mm[valid_linear] - error_threshold,
                     datum$obsYear[valid_linear], datum$linear_meas_mm[valid_linear] + error_threshold,
                     length = 0.05, angle = 90, code = 3, col = "gray80", lwd = 1.5)
    points(datum$obsYear, datum$linear_meas_mm, col = status[datum$QC + 2], pch = 17, cex = 2)
  }
  
  text(ulhc_x, ulhc_y, "FAIL\nALL", col = "red", font = 2, adj = c(0, 1))
  legend("bottomright", legend = c("fail", "unchecked", "pass", "Dendro Band (●)", "Linear Meas (▲)", "+/- 10mm Threshold"), 
         col = c(status[1:3], "black", "black", "gray80"), pch = c(18, 18, 18, 20, 17, NA), 
         lty = c(NA, NA, NA, NA, NA, 1), lwd = c(NA, NA, NA, NA, NA, 2), cex = 1.0, bty = "n")
  
  flag <- graphics::identify(click_x, click_y) 
  
  if (length(flag) > 0) {
    if (max(flag) > nrow(id_map)) {
      datum$QC[] <- -1
    } else {
      selected_rows <- id_map$row_idx[flag]
      datum$QC[selected_rows] <- -1
    }
  }
  
  # Transition remaining untouched values to passed (status 1)
  datum$QC[datum$QC == 0] <- 1
  
  # ==========================================================
  # PLOT 2: UNDO / REVIEW CAPABILITY
  # ==========================================================
  
  cat("\n------------------------------------------------------------\n")
  cat("QC CONSOLE INSTRUCTIONS - STEP 2 (UNDO/REVIEW):\n")
  cat("1. Review your choices. Outliers are marked in RED, accepted steps are BLUE.\n")
  cat("2. Click LEFT on any red outlier point to UNDO the flag (reverts it to 'pass').\n")
  cat("3. When satisfied, press ESC or RIGHT-CLICK and select 'Stop' to proceed to the next band.\n")
  cat("------------------------------------------------------------\n\n")
  
  plot(NULL, NULL, 
       xlim = x_limits, ylim = y_limits,
       xlab = "Observation Year", ylab = "Increment Measurement (mm)",
       main = paste("UPDATED:", meta_title, "\nClick points to UNDO outliers"))
  
  if (!all_dendro_na) {
    valid_dendro <- !is.na(datum$dendroRead_mm)
    graphics::arrows(datum$obsYear[valid_dendro], datum$dendroRead_mm[valid_dendro] - error_threshold,
                     datum$obsYear[valid_dendro], datum$dendroRead_mm[valid_dendro] + error_threshold,
                     length = 0.05, angle = 90, code = 3, col = "gray80", lwd = 1.5)
    points(datum$obsYear, datum$dendroRead_mm, col = status[datum$QC + 2], pch = 20, cex = 2)
  }
  if (!all_linear_na) {
    valid_linear <- !is.na(datum$linear_meas_mm)
    graphics::arrows(datum$obsYear[valid_linear], datum$linear_meas_mm[valid_linear] - error_threshold,
                     datum$obsYear[valid_linear], datum$linear_meas_mm[valid_linear] + error_threshold,
                     length = 0.05, angle = 90, code = 3, col = "gray80", lwd = 1.5)
    points(datum$obsYear, datum$linear_meas_mm, col = status[datum$QC + 2], pch = 17, cex = 2)
  }
  
  legend("bottomright", legend = c("fail", "unchecked", "pass"), 
         col = status[1:3], pch = 18, cex = 1.1, bty = "n")
  
  undo_flag <- graphics::identify(datum$obsYear[id_map$row_idx], id_map$value)
  
  if (length(undo_flag) > 0) {
    undo_rows <- id_map$row_idx[undo_flag]
    datum$QC[undo_rows] <- 1
  }
  
  return(invisible(datum))
}