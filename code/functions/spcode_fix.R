#' Interactive Manual Fix for Missing Species Taxonomy Codes
#'
#' This function scans a demographic summary dataset for missing taxonomy codes 
#' (\code{Sp_code}), surfaces the original field log records (\code{Sp_code.current}), 
#' and opens an interactive console interface allowing the researcher to manually 
#' input corrected species codes, genera, and specific epithets.
#'
#' @param data A data frame or tibble containing at least the following columns: 
#'   \code{Sp_code} (character), \code{Sp_code.current} (character), 
#'   \code{genus} (character), and \code{species} (character).
#'
#' @details The function loops sequentially through all row indices where 
#' \code{is.na(data$Sp_code)} evaluates to \code{TRUE}. For each row, it prints 
#' the unparsed field entry to the console and pauses execution using 
#' \code{\link{readline}} to await user input. 
#' 
#' If a user hits \code{Enter} without typing an entry, the field is safely 
#' preserved as \code{NA}.
#'
#' @return Returns the original data frame or tibble with updated strings in the 
#'   \code{Sp_code}, \code{genus}, and \code{species} columns for the modified rows.
#'
#' @author Your Name
#' 
#' @seealso \code{\link[stringr]{str_to_sentence}}
#' 
#' @examples
#' \dontrun{
#' # Assuming species_df has rows where Sp_code is NA
#' cleaned_species_df <- interactive_species_fix(species_df)
#' }
#'
#' @importFrom dplyr mutate
#' @export
Sp_code.fix <- function(data) {
  # Identify the indices of rows where Sp_code is missing
  na_indices <- which(is.na(data$Sp_code))
  
  if (length(na_indices) == 0) {
    message("No missing Sp_code values found! Data is clean.")
    return(data)
  }
  
  message(paste("Found", length(na_indices), "rows with missing Sp_code. Starting interactive entry...\n"))
  
  # Loop through each row that has an NA
  for (i in na_indices) {
    cat("--------------------------------------------------\n")
    cat("Row index:", i, "\n")
    cat("Field Notes Code (Sp_code.current):", data$Sp_code.current[i], "\n\n")
    
    # 1. Prompt for new Sp_code
    new_sp_code <- readline(prompt = "Enter new Sp_code (e.g., 'Acr acu'): ")
    
    # 2. Prompt for Genus
    new_genus <- readline(prompt = "Enter Genus name (e.g., 'Acrocomia'): ")
    
    # 3. Prompt for Species
    new_species <- readline(prompt = "Enter Species name (e.g., 'aculeata'): ")
    
    # Assign the inputs back to the data frame
    data$Sp_code[i] <- ifelse(new_sp_code == "", NA, new_sp_code)
    data$genus[i]   <- ifelse(new_genus == "", NA, new_genus)
    data$species[i] <- ifelse(new_species == "", NA, new_species)
    
    cat("\nSaved! Moving to next record...\n")
  }
  
  cat("--------------------------------------------------\n")
  message("Interactive entry complete!")
  return(data)
}