# Script to clean the flow cytometry data (remove negative controls, calculate cell counts, etc)
# Kendra E Walters
# May 12th, 2021

# Room of ()ments
require(data.table)
require(tidyverse)
require(stringr)

# Prepare data locations
# input.data <- read.csv("Kendra_processed_data/flow_cytometer_output__adjusted_gates.csv") %>% as_tibble()
## updated on 01/31/2025 to change how Kendra dealt with negatives. Instead of subtracting the average cell # from all samples, subtract the flow negatives from all samples and dispersal negatives (the ones out in the field) from the death slides only
## the flow negatives from the day of measurement were already subtracted from the samples in excel, these values are in the abs_count_minus_negative column
input.data <- read.csv("Kendra_processed_data/flow_cytometer_output__adjusted_gates_KB.csv") %>% as_tibble()
input.metadata <- load("Kendra_processed_data/metadata_cleaned.rda")

# Setting up here the samples that we want to keep (everything else is either a duplicate or was a negative control to make sure the flow cytometer is working but not something that we need to deal with in the analysis stage)
sample.keep <- c("n1", "n2", "n3", "n4", "g1 1:10", "g2 1:10", "g3 1:10", 
                 "g4 1:10", "g5 1:10", "g6 1:10",  "s1",  "s2 1:2", "s3 1:2", 
                 "s4 1:2", "s5 1:2","s6 1:2", "36 1:10(1)", "35", "34 1:10(2)", 
                 "33", "32 1:2", "31 1:2", "30 1:2", "29 1:5(1)", "28 1:10", 
                 "27 1:10", "26", "25 1:2", "24 1:20", "23 1:10", "22 1:10", 
                 "21 1:5(1)", "20 1:20", "19 1:50(1)", "18 1:10(2)", 
                 "17 1:10", "16", "15 1:5", "14 1:10(1)", "13 1:5", "12 1:5", 
                 "11 1:10", "10 1:10", "9 1:5", "8 1:10", "7 1:2", "6 1:5", 
                 "5", "4 1:5", "3 1:10", "2", "1 1:10(1)", "37 1:10", "38", 
                 "39", "40", "41 1:20", "42 1:10", "43 1:5", "45 1:10", "46 1:10", 
                 "47 1:10", "48 1:10", "49 1:10", "50  1:10", "51 1:20(1)", 
                 "52 1:50", "53 1:10", "54 1:5", "55 1:10", "56", "57 1:10", 
                 "58 1:20", "59 1:10", "60 1:20 from 1x", "61 1:10", "62 1:10", 
                 "63 1:20", "64 1:2", "65 1:10", "66 1:10 from 1x", "67 1:10", 
                 "68 1:10", "69 1:10 from 1x", "70 1:10", "71 1:25 from 1x", 
                 "72 1:5", "73 1:10", "74 1:10", "75 1:20", "76 1:20", 
                 "77 1:20", "78 1:10", "79 1:50", "80 1:5", "81 1:50", 
                 "82 1:10", "83 1:20", "84", "85 1:10 from OG", "87 1:20", "88 1:10", 
                 "89 1:10", "90 1:5", "91 1:50", "92 1:20", "93(1)", "94 1:20", 
                 "95 1:10 from OG", "96", "97 1:5", "98 1:10", "99 1:20", 
                 "100 1:20", "102", "103 1:20", "104", "105 1:10", "106 1:5", 
                 "107 1:10", "108 1:10", "109 1:10", "110 1:5", "111 1:10", 
                 "86 1:10", "112 1:10", "113", "114", "115 1:10", "116 1:10", 
                 "117 1:5", "118 1:10", "119 1:10", "120 1:10", "121 1:2", 
                 "122  1:10", "123 1:10", "124 1:20", "125 1:10", "126 1:10", 
                 "127 1:10", "128 1:20", "129 1:10", "130 1:5", "131", "132 1:10", 
                 "133 1:10", "134 1:10", "135 1:10", "136 1:10", "137 1:10", "138 1:10", 
                 "139 1:10  from OG", "140 1:10", "141 1:10", "142", "143 1:10",
                 "144", "145 1:2", "146 1:10", "147", "148 1:10", "149")


# Now to load the data and turn flow cytometry data into # of cells / glass slide and get our sample ID code
cell.counts <- input.data %>% # Load data
  filter(Sample %in% sample.keep) %>% # keep just those samples that I IDed on the day of running samples on the flow cytometer
  mutate(Dilution = str_match(Sample, ":([0-9]*)")[,2]) %>% # get dilution ratio used for each sample
  mutate(Dilution = as.numeric(Dilution)) %>% 
  mutate_if(is.numeric, replace_na, replace = 1) %>% 
  mutate(Starting_Volume_uL = 2000) %>% # the volume of saline added to each sample
  mutate(Mid_Volume_uL = ifelse(Specimen == "T0__Death__61920" &
                                  Sample %in% c("n1", "n2", "n3", "n4"), 
                                2000, 1500)) %>% # records any aliquots removed from saline + sample 
  mutate(End_Volume_uL = ifelse(Specimen == "T0__Death__61920" &
                                  Sample %in% c("n1", "n2", "n3", "n4"), 
                                2222, 1667)) %>% # records volume of 10% Pi-buffered GTA added to sample
  mutate(number_cells_on_glass_slide = Abs_count_minus_negative * Dilution * End_Volume_uL * (1/Mid_Volume_uL) * Starting_Volume_uL) %>% # uses these volumes to calculate the number of cells on the glass slide (before accounting for our negative controls)
  mutate(`Random_#_Flow` = str_extract(Sample, "[ngs]?[0-9]*"))


## Now we need to use our negative controls (laboratory negatives = n1, n2, ...; and field negatives = LC samples) to subtract any contamination that might have wormed its way into our samples
## On 01/31/2025 - KB changed this to only subtract the avg abundance on the negative samples in the field from the death slides. The immigration rate slides (open) do not need this done to them

# First step here -- calculate the mean contamination for the lab (just for T0 samples because they don't have any field negatives since they happened before the field!)
T0.neg.mean <- cell.counts %>% 
  filter(str_detect(Sample, "^n") & str_detect(Specimen.ID, "23")) %>% 
  summarize(mean(number_cells_on_glass_slide)) %>% as.numeric()

# Now calculating the mean contamination for the field controls (which would also include lab contamination so we don't need to subtract that separately)
metadata <- metadata.together %>% 
  select(c(New_Sample_ID, `Random_#_Flow_Cytometry`)) %>% # we just need these two columns today
  mutate(`Random_#_Flow`= as.character(`Random_#_Flow_Cytometry`)) # and changing the form here so it will join with our data

all.neg.means <- cell.counts %>% 
  left_join(metadata) %>% # joining the metadata and cell.counts dataframes (tibbles, I guess?? I can't take that name seriously)
  filter(str_detect(New_Sample_ID, "^LC")) %>%  # now we just want the glass slide samples that are Closed to immigration (negative controls)
  group_by(Specimen) %>% # group by timepoint
  summarize(mean_neg = mean(number_cells_on_glass_slide)) %>% # and then calculate the mean!
  bind_rows(tibble(Specimen = "T0__Death__61920", mean_neg = T0.neg.mean))

## added by Kristin on 01/30/2024
## calculate the overall mean for the in field negative controls per cm^2
closed.mean <- cell.counts %>% 
  left_join(metadata) %>% # joining the metadata and cell.counts dataframes (tibbles, I guess?? I can't take that name seriously)
  filter(str_detect(New_Sample_ID, "^LC")) %>%  # now we just want the glass slide samples that are Closed to immigration (negative controls)
  mutate(mean_cells_per_cm2 = (number_cells_on_glass_slide / (7.5 * 2.5))) %>% # convert counts from per glass slide to per cm^2
  summarize(mean_neg = mean(mean_cells_per_cm2), ## calculate mean
            sd = sd(mean_cells_per_cm2),
            n = sum(!is.na(mean_cells_per_cm2)),
            se = sd/sqrt(n)) ## calculate standard error

## calculate the mean abundance on open bags in the field
open.mean <- cell.counts %>% 
  left_join(metadata) %>% # joining the metadata and cell.counts dataframes (tibbles, I guess?? I can't take that name seriously)
  filter(str_detect(New_Sample_ID, "^LO")) %>%  # now we just want the glass slide samples that are Closed to immigration (negative controls)
  mutate(mean_cells_per_cm2 = (number_cells_on_glass_slide / (7.5 * 2.5))) %>% # convert counts from per glass slide to per cm^2
  summarize(mean_neg = mean(mean_cells_per_cm2), ## calculate mean
            sd = sd(mean_cells_per_cm2),
            n = sum(!is.na(mean_cells_per_cm2)),
            se = sd/sqrt(n)) ## calculate standard error

# ## back to kendra's script
# # And now the fun part! Subtracting these values (one mean for each T0 - T5) from our samples
# cleaned.cell.counts <- cell.counts %>% 
#   left_join(metadata) %>%
#   left_join(all.neg.means) %>% 
#   mutate(number_cells_on_glass_slide_cleaned = number_cells_on_glass_slide - mean_neg) %>% 
#   filter(str_detect(Sample, "^n", negate = TRUE))
# ## And write out our files!
# write_csv(cleaned.cell.counts, output.data)
# save(cleaned.cell.counts, file = output.rda)

## KB update - now just subtract the avg. negative control from the death rate samples only
cleaned.cell.counts.KB <- cell.counts %>% 
  left_join(metadata) %>% 
  left_join(all.neg.means) %>% 
  mutate(number_cells_on_glass_slide_cleaned = ifelse(grepl("LD", New_Sample_ID), 
                                                      number_cells_on_glass_slide - mean_neg,
                                                      number_cells_on_glass_slide))
write_csv(cleaned.cell.counts.KB, "cleaned_flow_data_KB.csv")

