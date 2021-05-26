# Analyze the contamination SourceTracker for 16S and ITS
# Kendra E Walters
# May 26th 2021


# Room of ()ments
require(tidyverse)
require(data.table)

# Load our cleaned data so we know what samples are over the rarefied level
load("03_Processed_data/16S_glass_slides.rda")
load("03_Processed_data/ITS_glass_slides.rda")

# Load mixing properties and label whether the samples were kept
contam.16S <- fread("01_Raw_data/03_Microbial_Community_Composition/glass_slides__16S/source_tracker_contamination/mixing_proportions.txt", data.table = FALSE) %>% 
  rename(Neg_proportion = Neg, 
         Unknown_proportion = Unknown) %>% 
  left_join(fread("01_Raw_data/03_Microbial_Community_Composition/glass_slides__16S/source_tracker_contamination/mixing_proportions_stds.txt", data.table = FALSE)) %>% 
  rename(Neg_std = Neg, 
         Unknown_std = Unknown) %>% 
  mutate(Kept = ifelse(SampleID %in% row.names(q1.rarefied), "Yes", "No"))

contam.ITS <- fread("01_Raw_data/03_Microbial_Community_Composition/glass_slides__ITS/source_tracker_contamination/source_tracker_contamination/mixing_proportions.txt", data.table = FALSE) %>% 
  rename(Neg_proportion = N, 
         Unknown_proportion = Unknown) %>% 
  left_join(fread("01_Raw_data/03_Microbial_Community_Composition/glass_slides__ITS/source_tracker_contamination/source_tracker_contamination/mixing_proportions_stds.txt", data.table = FALSE)) %>% 
  rename(Neg_std = N, 
         Unknown_std = Unknown) %>% 
  mutate(Kept = ifelse(SampleID %in% row.names(q1.ITS.rarefied), "Yes", "No"))

# Thoughts
# ITS = much lower contamination
# 16S = everything looks okay, despite higher proportions, because most samples with high proportions were taken out
# since we don't know where the contamination came from (could have been from the sequencer), we can just move forward here, knowing that the vast majority of samples are not matching with our negative controls