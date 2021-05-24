# Combine metadata files (I don't know why I have more than one but it's a problem)
# Kendra E Walters
# May 12 2021

# Room of ()ments
require(tidyverse)
require(data.table)

# Set up where we are pulling data from
setwd("01_Raw_data/01_metadata_files/")

input.data <- "01_Primary_Metadata - Sheet1.csv"
out.data <- "../../03_Processed_data/metadata_cleaned.csv"
out.rda <- "../../03_Processed_data/metadata_cleaned.rda"

# Read in the metadata for each timepoint, bind together, and clean up a bit
files <- list.files(".", pattern = "02_T.*csv")
timepoints <- map_dfr(files, fread) %>%  # read in and bind together
  select(-c(`Collected?`, V14)) %>% # clean up unneeded columns
  rename(`Random_#_FLow_Cytometry` = `Random_#`) # rename column to be more specific

# Read in the other metadata that has all timepoint samples, merge all together, and then clean up a bit
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
} # this lovely solution to my problem was taken from https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows/34096575#34096575

metadata.together <- fread(input.data) %>% # reading in data
  full_join(timepoints) %>% # merging the two metadata files
  mutate_cond(grepl("^[SG]_", New_Sample_ID), 
              Treatment = "Initial Litter",
              Transect = NA, 
              Location = NA, 
              Ecosystem = NA) %>% # filling in metadata for initial litter
  mutate_cond(grepl("Soil", New_Sample_ID), 
              Treatment = "Soil") %>% # filling in metadata for soil samples
  mutate_cond(grepl("^Env[GS]", New_Sample_ID), 
              Treatment = "Environmental Litter") %>% # filling in metadata for environmental litter
  mutate(Notes = ifelse(`Random_#_Extractions_Tubes` == "135", "Is actually two samples (135 + 127) combined", Notes)) %>% 
  mutate(Notes = ifelse(`Random_#_Extractions_Tubes` == "21", "Might be sample 21 + the positive PCR control?", Notes)) %>% 
  select(!c(Transect_Loc ,Transect_Loc...8, Transect_Loc...2, `Random_#_FLow_Cytometry`)) %>% 
  mutate(Transect_Loc = paste0(Transect, Location)) %>% 
  filter(!(is.na(Old_Sample_ID))) %>% 
  rename(`Random_#_Flow_Cytometry` = `Random_#_Flow`)

# to fix
# fix duplicated samples, pay attention to which one has the right timepoint (which may NOT match sampleID)


  #filter(!(New_Sample_ID == "LOT3A4" & is.na(Old_Sample_ID))) %>% View


# Write out the ONE metadata table (*chorus chimes*)
write_csv(metadata.together, out.data)
save(metadata.together, file = out.rda)  
  