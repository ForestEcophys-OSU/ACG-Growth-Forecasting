#'
#' # Dry Forest Growth Project
#' 
#' ## Script 1: Data cleaning
#'
#' ### Author: German Vargas G.
#' 
#' ### Date: 2026-05-18
#' 
# Load packages ------
library(tidyverse)

# Load functions ----
source("code/functions/dendro_check.R")# Function to check dendrometer readings
source("code/functions/spcode_fix.R")# Function to check Sp_code

# Load data -----
tree.dbh <- read.csv(file = "data/20260514_dendrometer_data_flagged.csv")

# 1. Replace NA in QC with 0 (not checked)
tree.dbh <- tree.dbh %>% 
  mutate(QC = ifelse(test = is.na(QC),yes = 0,no = QC))
head(tree.dbh)

# Quality check dendrometer readings -----
#data.per.plot <- list()
load(file = "data/dendrocheckList.R")
lapply(data.per.plot, head)

data.per.plot[[21]] <- tree.dbh %>% 
  filter(plot == 21) %>%
  arrange(plot,subplot,dendrometer,dendroNumber,obsYear) %>%
  group_by(plot,subplot,dendrometer,dendroNumber) %>%
  do(dendro.check(.)) %>%
  ungroup() %>% 
  as.data.frame()

save(data.per.plot,file = "data/dendrocheckList.R")

tree.dbh2 <- do.call("rbind",data.per.plot) %>%
  select(plot,subplot,tree,stem,dendrometer,Sp_code,obsYear,DBH_tape_cm,
         dendroRead_mm,linear_meas_mm,liana,phenology,canopy,Notes,mortality,
         newDendrometer,dendroNumber,QC) %>%
  mutate(linear_meas_mm = ifelse(linear_meas_mm == 0.0,yes = NA,no = linear_meas_mm),# Set any linear measurements that are 0 to NA for cleaner formatting
         Notes = ifelse(Notes %in% "", yes = NA,no = Notes),# Set any empty Notes values to NA for cleaner formatting
         QC = ifelse(Notes %in% "not found", yes = -1,no = QC),# If "not found" is found anywhere in the Notes string, QC becomes -1
         QC = if_else(condition = !is.na(Notes) & str_detect(Notes, "dendro read error"), true = -1, false = QC),# If "dendro read error" is found anywhere in the Notes string, QC becomes -1
         Notes = if_else(condition = is.na(Notes) & QC == -1, true = "dendro read error", false = Notes)) %>%# If QC is -1 but Notes is NA, we will set Notes to "dendro read error" for easier identification of records with dendrometer reading errors that need further review.
  as.data.frame()

# Check species codes -------
sp_codes <- read.csv(file = "data/20211201_SpeciesList.csv")

# 1 - Fix typos in names first 
tree.dbh2 <- tree.dbh2 %>%
  mutate(Sp_code = str_replace(Sp_code, "^([A-Za-z]{3})([A-Za-z]{3})$", "\\1 \\2"), # Add space between the first three and last three letters in the species code, which is the format used in our species list#
         Sp_code = str_to_sentence(Sp_code), # Convert species codes to sentence case to match the format in our species list
         Sp_code = if_else(condition = Sp_code %in% "Dio spi", true = "Dio sal", false = Sp_code), # fix code for Dio sal, which was misspelled as Dio spi in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Dis ame", true = "Dip ame", false = Sp_code), # fix code for Dip ame, which was misspelled as Dis ame in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Eug man", true = "Eug mon", false = Sp_code), # fix code for Eug mon, which was misspelled as Eug man in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Fic sp", true = "Fic spp", false = Sp_code), # fix code for Fic spp, which was misspelled as Fic sp in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Paq qui", true = "Bom qui", false = Sp_code), # fix code for Bom qui, which was misspelled as Pac qui in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Pit duc", true = "Pit dul", false = Sp_code), # fix code for Pit dul, which was misspelled as Pit duc in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Tev ova", true = "The ova", false = Sp_code), # fix code for The ova, which was misspelled as Tev ova in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Tri syl", true = "Tri spp", false = Sp_code), # fix code for The spp, which was misspelled as Tri syl in the original dataset
         Sp_code = if_else(condition = Sp_code %in% "Vac col", true = "Aca col", false = Sp_code), # fix code for Aca col, which was classified as Vac col in the original dataset - will handle update on species name later
         Sp_code = if_else(condition = Sp_code %in% "", true = NA, false = Sp_code)) # Set any empty species codes to NA for further review

current_sp_codes <- tree.dbh2 %>%
  group_by(Sp_code) %>% 
  summarise(n = n()) %>% # Check all species codes in the dataset, and how many records they have. This will help us identify any potential issues with species codes, such as typos or inconsistencies.
  mutate(Sp_code.current = Sp_code,# Create a new column to store the current species code for comparison
                            Sp_code = if_else(condition = Sp_code.current %in% sp_codes$Sp_code, # Check if the current species code is in the list of valid species codes from our species list
                                              true = Sp_code.current, # If it is valid, keep it as is 
                                              false = NA)) %>% # otherwise, set it to NA for further review
  left_join(sp_codes %>% select(Sp_code,genus,species), by = c("Sp_code" = "Sp_code")) %>% # Join with the species list to get the corresponding genus and species names
  print(n = 130)# Print the updated species codes to check which ones were flagged as NA for review

current_sp_codes.clean <- Sp_code.fix(current_sp_codes) # Use the Sp_code_fix function to evaluate missing Sp_code in the species list data record
current_sp_codes.clean %>% print(n = 130) # Print the cleaned species codes to check which ones were corrected and which ones still need review
rm(current_sp_codes, sp_codes)

# Updated taxonomic nomenclature in the WFO backbone data ----
#' For this we will use the World Flora Online (WFO) backbone data, which is a comprehensive and authoritative source of plant taxonomic information. 
#' The WFO provides a standardized list of plant names and their corresponding taxonomic details, which can be used to verify and correct species codes 
#' in our dataset.
#' 
#' Citation: The World Flora Online Consortium, Elliott, A., Hyam, R., Watson, M., Wrankmore, E., Hartley, H., Krieger, J., Gandhi, K., 
#' Abad-Brotons, J., Acuña, R., Alcantara, S., Almeida, R. F. D., Alonso-Vargas, M. Á., Amorim, G., Anderson, G., Andrella, G. C., 
#' Anguiano, M., Antonio-Domingues, H., Ardi, W. H., … Španiel, S. (2025). World Flora Online Plant List December 2025 (2025-12) [Data set]. 
#' Zenodo. https://doi.org/10.5281/zenodo.18007552
#' 
library(WorldFlora)
#WFO.download(save.dir = "data/WFO_dir", version = "2025-12") # no need to run until next update, as we have the latest version downloaded
WFO.remember(WFO.file = "data/WFO_dir/classification.csv") # Load the WFO backbone data into R for use in taxonomic verification and correction

current_sp_codes.clean <- current_sp_codes.clean %>%
   mutate(spec.name.ORIG = if_else(condition = is.na(species), true = genus, false = paste(genus, species, sep = " "))) # Create a new column that combines the genus and species names for easier comparison with the WFO backbone data
current_sp_codes.clean %>% print(n = 130) # Print the updated species codes with the new GenusSpecies column to check for any remaining issues
colnames(current_sp_codes.clean) <- c("Sp_code.ORIG","n","Sp_code.current","genus.ORIG","species.ORIG","spec.name.ORIG")# change names to match the format used in the WFO.match function

# Use the WFO.match function to compare the spec.name.ORIG names in our dataset with the WFO backbone data, and identify any discrepancies or 
# updates needed for taxonomic nomenclature. This will help us ensure that our species codes are accurate and up-to-date according to the latest taxonomic standards.
wfo_report <- WFO.match(spec.data = current_sp_codes.clean$spec.name.ORIG,WFO.data = WFO.data,acceptedNameUsageID.match = T)

# Review the WFO.match report to identify any discrepancies or updates needed for taxonomic nomenclature. 
wfo_report_clean <- wfo_report %>% 
  distinct(spec.name.ORIG, spec.name, family, genus, specificEpithet, scientificName, .keep_all = TRUE) %>%
  select(spec.name.ORIG, spec.name, family, genus, specificEpithet, scientificName,scientificNameAuthorship,Old.status,taxonomicStatus,New.accepted) %>%
  filter(!scientificName %in% c("Annona mucosa","Byrsonima verbascifolia","Cecropia pachystachya","Dalbergia pervillei",
                                "Hirtella triandra subsp. triandra","Tabebuia ochracea","Premna thwaitesii","Premna micrantha"))

# Join the cleaned species codes with the WFO.match report to create a final species list that includes the updated taxonomic nomenclature and any necessary corrections to the species codes. 
final_sp_list <- current_sp_codes.clean %>% # Start with the cleaned species codes from our dataset
  left_join(wfo_report_clean, by = c("spec.name.ORIG" = "spec.name.ORIG")) %>% # Join with the cleaned WFO.match report to get the updated taxonomic information for each species
  select(Sp_code.current, genus.ORIG, species.ORIG, spec.name.ORIG, # Select relevant columns for the final species list
         family, genus, specificEpithet, scientificName, scientificNameAuthorship, 
         Old.status, taxonomicStatus, New.accepted) %>% 
  rename(Sp_code.ORIG = Sp_code.current) %>% # Rename the current species code column to Sp_code.ORIG for clarity in comparison with the updated species code we will generate based on the WFO backbone dat
  mutate(specificEpithet = if_else(condition = specificEpithet == "", true = NA, false = specificEpithet),# Set any empty specific epithet values to NA for cleaner formatting
         # Set any empty Old.status values to NA for cleaner formatting
         Old.status = if_else(condition = Old.status == "", true = NA, false = Old.status), 
         # If the scientific name is NA, we will set the family, genus, specific epithet, and scientific name to "unknown" for easier identification of records that need further review.
         across(.cols = c(family, genus, specificEpithet, scientificName),
                .fns = ~ if_else(condition = is.na(scientificName), true = "unknown", false = .x)),
         # If the specific epithet is unknown, we cannot generate a species code based on the genus and specific epithet, so we will set the species code to "unknown" as well. 
         # Otherwise, we will generate a new species code based on the first three letters of the genus and specific epithet, which is the format used in our species list. 
         Sp_code = if_else(condition = specificEpithet == "unknown",
                           true = "unknown", false = paste(str_sub(string = genus,start = 1,end = 3), str_sub(string = specificEpithet,start = 1,end = 3)),
                           missing = paste(str_sub(string = genus,start = 1,end = 3), "spp")),
         # If the current species code in our dataset is NA, we will set it to "unknown" for easier identification of records that need further review. 
         # Otherwise, we will keep the current species code as is for comparison with the updated species code based on the WFO backbone data.
         Sp_code.ORIG = if_else(condition = is.na(Sp_code.ORIG),true = "unknown",false = Sp_code.ORIG)) %>%
  arrange(family, genus, specificEpithet) %>%
  print(n = 130)

# Save the final species list with the updated taxonomic nomenclature and any necessary corrections to the species codes for future reference and use in our dataset.
write_csv(final_sp_list, file = "data/20260526_SpeciesList.csv", na = "NA") 

# Update the species codes in our dataset based on the final species list -----
tree.dbh2 <- tree.dbh2 %>% 
  rename(Sp_code.ORIG = Sp_code) %>%
  mutate(Sp_code.ORIG = if_else(condition = is.na(Sp_code.ORIG),true = "unknown",false = Sp_code.ORIG),
         Notes = if_else(condition = is.na(subplot), true = "check subplot",false = Notes),
         Notes = if_else(condition = Sp_code.ORIG == "unknown", true = "check sp", false = Notes)) %>%
  left_join(final_sp_list %>% select(Sp_code.ORIG, Sp_code), by = c("Sp_code.ORIG" = "Sp_code.ORIG")) %>%
  select(plot,subplot,tree,stem,dendrometer,Sp_code,Sp_code.ORIG,obsYear,DBH_tape_cm,
         dendroRead_mm,linear_meas_mm,liana,phenology,canopy,Notes,mortality,
         newDendrometer,dendroNumber,QC)

# Save latest version of data -----
write_csv(tree.dbh2,file = "data/20260626_dendrometer_data_flagged.csv",na = "")
