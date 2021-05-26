# Script to make sure our sequencing looks good, then remove any controls, mock communities, etc
# Output: rarefied and rounded OTU table
# Output: alpha-diversity (Shannon, richness, and Simpson) metrics for each sample in a single table
# Output: average Bray-Curtis matrix rarefied to sample depth as rarefied OTU table
# 16S amplicon sequencing on glass slides
# Kendra E Walters
# May 20th 2021

## Room of ()ments
require(data.table)
require(tidyverse)
require(vegan)
require(EcolUtils)


## Setting up inputs and outputs
setwd("01_Raw_data/03_Microbial_Community_Composition/glass_slides__16S/")

source.tracker.file <- "16S_glass_slides__for_source_tracker.tsv"
source.tracker.meta <- "map_16S_glass_slides__sourcetracker.txt"
source.tracker.q1.file <- "../../../03_Processed_data/16S_glass_slides_q1__source_tracker.tsv"
source.tracker.q1.meta <- "../../../03_Processed_data/16S_glass_slides_q1__source_tracker_map.tsv"

output.rarefied.table <- "../../../03_Processed_data/16S_glass_slides__rarefied_971_rounded.tsv"
output.bc.matrix <- "../../../03_Processed_data/16S_glass_slides__BC_median__r971.tsv"
output.rda <- "../../../03_Processed_data/16S_glass_slides.rda"

## Setting up parameters
rarefaction.depth <- 971


## Load data
q1.table <- as.data.frame(fread("12_table_filtered.tsv")) %>% as_tibble() %>% 
  select(!(c(LDT2R10CG, LOT4A1))) # these two samples were contaminated in the PCR step so I'm removing them here
#   --p-trim-left 5 \
#   --p-trunc-len 224 \
q1.combined <- as.data.frame(fread("14_table_filtered_with_taxonomy_filtered.tsv")) %>% as_tibble() %>% 
  select(!(c(LDT2R10CG, LOT4A1))) # these two samples were contaminated in the PCR step so I'm removing them here

load("../../../03_Processed_data/metadata_cleaned.rda")


## First, we check the PCR+ controls (our mock communities)
q1.combined %>% select(contains("PCR+"), taxonomy, `#OTU ID`) %>% rowwise %>% 
  filter(sum(c_across(starts_with("PCR+"))) > 0) %>% ungroup %>% View
# is PCR+_9_16_20 actually a negative? It seems like it's not a positive control. Using it as a negative moving forward
# For the other samples, we see: Pseudomonas, Lactobacillus, Enterococcus, Listeriaceae, Bacillus, Salmonella, Staphylococcus, Salmonella enterica


## Next, let's check our PCR- controls 
q1.combined %>% select(contains("PCR-"), `PCR+_9_16_2`, taxonomy, `#OTU ID`) %>% rowwise %>% 
  filter(sum(c_across(starts_with("PCR"))) > 0) %>% ungroup %>% View
# These look pretty good, I would say! Really low sequence counts for all so they would just be filtered out by rarefaction if we kept them in


## Now let's check the extraction negatives (just reagents, no samples). These are super important for our low biomass samples
q1.combined %>% select(contains("Neg"), taxonomy, `#OTU ID`) %>% rowwise %>% 
  filter(sum(c_across(starts_with("Neg"))) > 0) %>% ungroup %T>% View %>% select(starts_with("Neg")) %>% colSums
# all of these are beneath our rarefaction level


## Setting up our feature table and metadata to use SourceTracker to check out contamination sources
data.frame(c("# Constructed from biom file")) %>% write_tsv(source.tracker.file, col_names = FALSE)
q1.source.tracker <- q1.table %>% select(!(contains("PCR"))) %>% write_tsv(source.tracker.file, append = TRUE, col_names = TRUE)

data.frame("SampleID" = names(q1.source.tracker %>% select(!(contains("OTU"))))) %>% as_tibble %>% 
  mutate(SourceSink = ifelse(grepl("Neg", SampleID), "source", "sink")) %>% 
  mutate(Env = str_match(SampleID,"^([a-zA-SU-Z]+)_?T?[0-9]")[,2]) %>% 
  rename("#SampleID" = SampleID) %>% 
  write_tsv(source.tracker.meta)

# Run the following lines in the command line
# biom convert -i 16S_glass_slides__for_source_tracker.tsv -o 16S_glass_slides__for_source_tracker.biom --to-hdf5 --table-type="OTU table"
# source activate st2
# sourcetracker2 gibbs -i 16S_glass_slides__for_source_tracker.biom -m map_16S_glass_slides__sourcetracker.txt --source_rarefaction_depth 0 --sink_rarefaction_depth 0  -o source_tracker_contamination/
# Note: this took over 48 hours to run on a single core on my lil MacBook




## Output: rarefied & rounded OTU table
# Decide on rarefaction depth
sums <- colSums(q1.table[, 2:ncol(q1.table)])
hist(sums[sums < 20000], breaks = 100)
hist(sums[sums < 5000], breaks = 100)

sort(sums) # maybe rarefy to 971?
length(sums[sums >=971]) / length(sums) # keep 63.8% of samples

# Make the OTU ID the row name, drop samples under rarefaction depth, and transform
# because vegan expects rows = samples and columns = species
q1.rarefied <- q1.table %>% 
  column_to_rownames('#OTU ID') %>% 
  select(which(colSums(.) >= rarefaction.depth)) %>% 
  select(!(contains("PCR"))) %>%
  select(!(contains("AirT4"))) %>% # these samples fell on the ground during the field experiment and aren't truly air samples
  t %>% as.data.frame() %>% 
  rrarefy.perm(., sample = rarefaction.depth, n = 1000, round.out = TRUE) %>% 
  as.data.frame

# Write file
#write.table(q1.rarefied, file = output.rarefied.table, quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)



## Output: alpha-diversity (Shannon, richness, and Simpson) metrics for each sample in a single table
q1.alpha <- data.frame("Shannon" = diversity(q1.rarefied, index = "shannon"), 
           "Richness" = apply(q1.rarefied[,]>0, 1, sum), 
           "Simpson" = diversity(q1.rarefied, index = "simpson"))
  


## Output: average Bray-Curtis matrix rarefied to sample depth as rarefied OTU table
q1.bray.dist <- q1.table %>% 
  column_to_rownames('#OTU ID') %>% 
  select(which(colSums(.) >= rarefaction.depth)) %>% # remove samples that fall under rarefaction depth
  select(!(contains("PCR"))) %>% # remove our PCR+ and - 
  select(!(contains("AirT4"))) %>% #these samples fell on the ground during the field experiment and aren't truly air
  t %>% as.data.frame() %>% 
  avgdist(., sample = 971, meanfun = median, transf = sqrt, iterations = 999) #default is Bray-Curtis calculated by vegdist
q1.bray.dist.df <- as.data.frame(as.matrix(q1.bray.dist))
  # write.table(x = ., output.bc.matrix, sep = "\t", row.names = FALSE)




## Check to make sure our PCR replicates (same samples run through PCR multiple times) overlap
dup.samples <- c("LDT0R5G", "LDT0R5G_oops") # these are low biomass samples so our other PCR replicates did not sequence
NMDS <- metaMDS(q1.bray.dist, autotransform = FALSE, k = 2, trymax = 200)
data.frame(NMDS$points[,1:2]) %>% 
  mutate(Duplicates = ifelse(row.names(.) %in% dup.samples, "Yes", "No")) %>% 
  ggplot(data = .) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Duplicates)), size = 4) + 
  theme_classic() +
  scale_color_brewer(palette = "Set1") 
# looks good!



## Save the pertinent objects as an R data file
save(metadata.together, q1.alpha, q1.rarefied, q1.bray.dist, q1.bray.dist.df, file = output.rda)



## Last, but not least, we will output the needed files to do a SourceTracker analysis to answer: 
## Where are the dispersing taxa coming from? (SourceTracker, dispersal community ~ air + soil + litter)
data.frame(c("# Constructed from biom file")) %>% write_tsv(source.tracker.q1.file, col_names = FALSE)
q1.rarefied %>% 
  filter(!(grepl("LD", row.names(.)))) %>% 
  rownames_to_column("SampleID") %>% as_tibble() %>%
  mutate(SampleID = ifelse(SampleID == "LOT1A3_dup_new_protocol", "LOT1A3", SampleID)) %>% 
  column_to_rownames("SampleID") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("#OTU ID") %>% 
  write_tsv(source.tracker.q1.file, append = TRUE, col_names = TRUE)

key <- c("1" = "Shrubland", "2" = "Shrubland", "3" = "Shrubland", "4" = "Shrubland",
         "5" = "Grassland", "6" = "Grassland", "7" = "Grassland", "8" = "Grassland")

data.frame("SampleID" = row.names(q1.rarefied %>% 
                                    filter(!(grepl("LD", row.names(.)))) %>% 
                                    rownames_to_column("SampleID") %>% as_tibble() %>%
                                    mutate(SampleID = ifelse(SampleID == "LOT1A3_dup_new_protocol", "LOT1A3", SampleID)) %>%
                                    column_to_rownames("SampleID"))) %>%  # just get rownames out 
  mutate(SourceSink = ifelse(grepl("Air|Soil|Env", SampleID), "source", "sink")) %>% 
  mutate(Env = str_match(SampleID,"^([a-zA-SU-Z]+)_?T?[0-9]")[,2]) %>% 
  mutate(Env = ifelse(grepl("LO", Env), 
                      key[str_match(SampleID, "[ABC]([1-8])")[,2]], # converting from place on transect to ecosystem
                      Env)) %>%
  rename("#SampleID" = SampleID) %>%
  write_tsv(source.tracker.q1.meta)
# Run:
# source activate st2
# biom convert -i 16S_glass_slides_q1__source_tracker.tsv -o 16S_glass_slides_q1__source_tracker.biom --to-hdf5 --table-type="OTU table"
# sourcetracker2 gibbs -i 16S_glass_slides_q1__source_tracker.biom -m 16S_glass_slides_q1__source_tracker_map.tsv --source_rarefaction_depth 0 --sink_rarefaction_depth 0 --jobs 2 -o 16S_glass_slide__dispersal_source_tracker/

