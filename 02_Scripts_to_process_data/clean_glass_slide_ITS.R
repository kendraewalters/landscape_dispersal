# Script to make sure our sequencing looks good, then remove any controls, mock communities, etc
# Output: rarefied and rounded OTU table
# Output: alpha-diversity (Shannon, richness, and Simpson) metrics for each sample in a single table
# Output: average Bray-Curtis matrix rarefied to sample depth as rarefied OTU table
# ITS amplicon sequencing on glass slides
# Kendra E Walters
# May 20th 2021

## Room of ()ments
require(data.table)
require(tidyverse)
require(vegan)
require(EcolUtils)


## Setting up inputs and outputs
setwd("01_Raw_data/03_Microbial_Community_Composition/glass_slides__ITS/")

source.tracker.contam.file <- "ITS_glass_slides_contamination__for_source_tracker.tsv"
source.tracker.contam.meta <- "map_ITS_glass_slides_contamination__sourcetracker.txt"
source.tracker.q1.ITS.file <- "../../../03_Processed_data/ITS_glass_slides_q1__source_tracker.tsv"
source.tracker.q1.ITS.meta <- "../../../03_Processed_data/ITS_glass_slides_q1__source_tracker_map.tsv"

output.rarefied.table <- "../../../03_Processed_data/ITS_glass_slides__rarefied_2400_rounded.tsv"
output.bc.matrix <- "../../../03_Processed_data/ITS_glass_slides__BC_median__r2400.tsv"
output.rda <- "../../../03_Processed_data/ITS_glass_slides.rda"

## Setting up parameters
rarefaction.depth <- 2400


## Load data
q1.ITS.table <- as.data.frame(fread("12_table_filtered.tsv")) %>% as_tibble()
#   --p-trim-left 5 \
#   --p-trunc-len 225 \
q1.ITS.combined <- as.data.frame(fread("14_table_filtered_with_taxonomy_filtered.tsv")) %>% as_tibble()

load("../../../03_Processed_data/metadata_cleaned.rda")


## First, we check the PCR+ controls (our mock communities)
q1.ITS.combined %>% select(contains("PCR+"), taxonomy, `#OTU ID`) %>% rowwise %>% 
  filter(sum(c_across(starts_with("PCR+"))) > 0) %>% ungroup %>% View
# Here, we see: Saccharomyces, Cryptococcus_neoformans, with multiple OTUs per ID
# All sequenced well except for PCR+_10_6_20, I guess that just happens sometimes


## Note: we did not have PCR- controls that got sequenced for ITS samples


## Now let's check the extraction negatives (just reagents, no samples). These are super important for our low biomass samples
q1.ITS.combined %>% select(N3, N4, N5, taxonomy, `#OTU ID`) %>% 
  rowwise %>% filter(sum(c_across(starts_with("N"))) > 0) %>% ungroup %T>% View %>% select(starts_with("N")) %>% colSums
# all of these are below 400 sequences, seems fairly decent


## Setting up our feature table and metadata to use SourceTracker to check out contamination sources
data.frame(c("# Constructed from biom file")) %>% write_tsv(source.tracker.contam.file, col_names = FALSE)
q1.ITS.source.tracker <- q1.ITS.table %>% select(!(contains("PCR"))) %>% write_tsv(source.tracker.contam.file, append = TRUE, col_names = TRUE)

data.frame("SampleID" = names(q1.ITS.source.tracker %>% select(!(contains("OTU"))))) %>% as_tibble %>% 
  mutate(SourceSink = ifelse(grepl("^N", SampleID), "source", "sink")) %>% 
  mutate(Env = str_match(SampleID,"^([a-zA-SU-Z]+)_?T?[0-9]")[,2]) %>% 
  rename("#SampleID" = SampleID) %>% 
  write_tsv(source.tracker.contam.meta)

# Run the following lines in the command line
# biom convert -i ITS_glass_slides_contamination__for_source_tracker.tsv -o ITS_glass_slides_contamination__for_source_tracker.biom --to-hdf5 --table-type="OTU table"
# source activate st2
# sourcetracker2 gibbs -i ITS_glass_slides_contamination__for_source_tracker.biom -m map_16S_glass_slides_contamination__sourcetracker.txt --source_rarefaction_depth 0 --sink_rarefaction_depth 0  -o source_tracker_contamination/
# Note: this took about 48 hours to run on a single core on my lil MacBook




## Output: rarefied & rounded OTU table
# Decide on rarefaction depth
sums <- colSums(q1.ITS.table[, 2:ncol(q1.ITS.table)])
hist(sums[sums < 20000], breaks = 100)
hist(sums[sums < 5000], breaks = 100)

sort(sums) # maybe rarefy to 971?
length(sums[sums >=2400]) / length(sums) # keep 80% of samples

# Make the OTU ID the row name, drop samples under rarefaction depth, and transform
# because vegan expects rows = samples and columns = species
q1.ITS.rarefied <- q1.ITS.table %>% 
  column_to_rownames('#OTU ID') %>% 
  select(which(colSums(.) >= rarefaction.depth)) %>% 
  select(!(contains("PCR") | starts_with("N"))) %>% 
  select(!(contains("AirT4"))) %>% # these samples fell on the ground during the field experiment and aren't truly air samples
  t %>% as.data.frame() %>% 
  rrarefy.perm(., sample = rarefaction.depth, n = 1000, round.out = TRUE) %>% 
  as.data.frame

# Write file
#write.table(q1.ITS.rarefied, file = output.rarefied.table, quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)



## Output: alpha-diversity (Shannon, richness, and Simpson) metrics for each sample in a single table
q1.ITS.alpha <- data.frame("Shannon" = diversity(q1.ITS.rarefied, index = "shannon"), 
                       "Richness" = apply(q1.ITS.rarefied[,]>0, 1, sum), 
                       "Simpson" = diversity(q1.ITS.rarefied, index = "simpson"))



## Output: average Bray-Curtis matrix rarefied to sample depth as rarefied OTU table
q1.ITS.bray.dist <- q1.ITS.table %>% 
  column_to_rownames('#OTU ID') %>% 
  select(which(colSums(.) >= rarefaction.depth)) %>% # remove samples that fall under rarefaction depth
  select(!(contains("PCR") | starts_with("N"))) %>% # remove our PCR+ and - 
  select(!(contains("AirT4"))) %>% # these samples fell on the ground during the field experiment and aren't truly air samples
  t %>% as.data.frame() %>% 
  avgdist(., sample = 971, meanfun = median, transf = sqrt, iterations = 999) #default is Bray-Curtis calculated by vegdist
q1.ITS.bray.dist.df <- as.data.frame(as.matrix(q1.ITS.bray.dist))
  #write.table(x = ., output.bc.matrix, sep = "\t", row.names = FALSE)




## Note: we did not run any PCR replicates because we ran out of DNA so each sample got run once, or not at all, with ITS PCR



## Save the pertinent objects as an R data file
save(metadata.together, q1.ITS.alpha, q1.ITS.rarefied, q1.ITS.bray.dist, q1.ITS.bray.dist.df, file = output.rda)



## Last, but not least, we will output the needed files to do a SourceTracker analysis to answer: 
## Where are the dispersing taxa coming from? (SourceTracker, dispersal community ~ air + soil + litter)
data.frame(c("# Constructed from biom file")) %>% write_tsv(source.tracker.q1.ITS.file, col_names = FALSE)
q1.ITS.rarefied %>% 
  filter(!(grepl("LD", row.names(.)))) %>% 
  t %>% as.data.frame() %>% 
  rownames_to_column("#OTU ID") %>% 
  write_tsv(source.tracker.q1.ITS.file, append = TRUE, col_names = TRUE)

key <- c("1" = "Shrubland", "2" = "Shrubland", "3" = "Shrubland", "4" = "Shrubland",
         "5" = "Grassland", "6" = "Grassland", "7" = "Grassland", "8" = "Grassland")

data.frame("SampleID" = row.names(q1.ITS.rarefied %>% 
                                    filter(!(grepl("LD", row.names(.)))) %>% 
                                    rownames_to_column("SampleID") %>% as_tibble() %>%
                                    column_to_rownames("SampleID"))) %>%  # just get rownames out 
  mutate(SourceSink = ifelse(grepl("Air|Soil|Env", SampleID), "source", "sink")) %>% 
  mutate(Env = str_match(SampleID,"^([a-zA-SU-Z]+)_?T?[0-9]")[,2]) %>% 
  mutate(Env = ifelse(grepl("LO", Env), 
                      key[str_match(SampleID, "[ABC]([1-8])")[,2]], # converting from place on transect to ecosystem
                      Env)) %>%
  rename("#SampleID" = SampleID) %>% 
  write_tsv(source.tracker.q1.ITS.meta)

# Run:
# source activate st2
# biom convert -i ITS_glass_slides_q1__source_tracker.tsv  -o ITS_glass_slides_q1__source_tracker.biom --to-hdf5 --table-type="OTU table"
# sourcetracker2 gibbs -i ITS_glass_slides_q1__source_tracker.biom -m ITS_glass_slides_q1__source_tracker_map.tsv --source_rarefaction_depth 0 --sink_rarefaction_depth 0 --jobs 2 -o ITS_glass_slide__dispersal_source_tracker/


