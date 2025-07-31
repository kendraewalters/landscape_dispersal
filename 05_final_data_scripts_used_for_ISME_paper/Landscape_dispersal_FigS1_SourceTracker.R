## Figure S1 SourceTracker Results in Grassland and Shrubland
## Created by KB on 12/12/2024 after re-processing all sequences through Qiime2 and last edited it on July 21st, 2025
## Updated from Kendra E Walters' script created on July 12 2021

library(vegan)
library(car)
library(data.table)
library(plyr)
library(dplyr)
library(devtools)
library(EcolUtils)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(stringr)
library(forcats)
library(RColorBrewer)

## Set working directory
setwd("~/kbarbou1@uci.edu - Google Drive/My Drive/Professional_development/Kendra_Ch3_Manuscript/03_Analysis/")

## import metadata
metadata.16S <- read.delim("KB_qiime2_redo/Qiime2_output/q1_16S_barcodes_and_pooling__metadata.tsv",  check.names = FALSE) 
metadata.ITS <- read.delim("KB_qiime2_redo/Qiime2_output/q1_ITS_barcodes_and_pooling__metadata.tsv", check.names = FALSE) %>% column_to_rownames(., var="#SampleID")

## remove controls (all negative controls and PCR mock communities look good ~ negatives are below rarefaction and mock communities match expected taxa)
## remove samples contaminated at the PCR step
metadata.16S.nc <- metadata.16S %>% filter(!Ecosystem=="NA") %>% column_to_rownames(., var="SampleID")
metadata.ITS.nc <- metadata.ITS %>% filter(!Ecosystem=="NA")

##### SourceTracker #####
## identify sources and sinks in the rarefied datasets
## first get a single rarefied table for 16S and ITS
## import cleaned OTU table
otu.table.clean.16S <- read.delim("otu.table.clean.16S.v1.tsv", sep="\t")
## rarefy
rared.otu.16S <- otu.table.clean.16S %>% 
  dplyr::select(which(colSums(.) >= 970)) %>% ## select samples that meet rarefaction cut off
  t() %>% 
  as.data.frame() %>% 
  rrarefy.perm(., sample = 970, n = 1000, round.out = T) %>% 
  as.data.frame()
sort(rowSums(rared.otu.16S))
## remove columns with less than 1 OTU
rared.otu.16S <- as.data.frame(rared.otu.16S[, colSums(rared.otu.16S) >= 1]) ## kept 114 samples
write.table(rared.otu.16S, file="rared_otu_16S.tsv", sep="\t")

## import cleaned OTU table
otu.table.clean.ITS <- read.delim("otu.table.clean.ITS.v1.tsv", sep="\t")
##rarefy
rared.otu.ITS <- otu.table.clean.ITS %>% 
  dplyr::select(which(colSums(.) >= 1840)) %>% ## select samples that meet rarefaction cut off
  t() %>% 
  as.data.frame() %>% 
  rrarefy.perm(., sample = 1840, n = 1000, round.out = T) %>% 
  as.data.frame()
sort(rowSums(rared.otu.ITS))
## keep only rows close to rarefaction level and remove columns with less than 1 OTU
rared.otu.ITS <- as.data.frame(rared.otu.ITS[, colSums(rared.otu.ITS) >= 1]) ## kept 120 samples
write.table(rared.otu.ITS, file="rared_otu_ITS.tsv", sep="\t")

## get files ready for sourcetracker, remove death slides and move rownames to separate column 
sourcetracker.OTU.16S <- rared.otu.16S %>% 
  .[!grepl("LD", rownames(.)), ] %>% 
  rownames_to_column(var="SampleID") %>% 
  mutate(SampleID = ifelse(SampleID == "LOT1A3_dup_new_protocol", "LOT1A3", SampleID)) %>% 
  column_to_rownames("SampleID") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("#OTU ID")
write.table(sourcetracker.OTU.16S, file = "sourcetracker.OTU.16S.tsv", sep="\t") ## add "# Constructed from biom file" in terminal

## get metadata files to match
key <- c("1" = "Shrubland", "2" = "Shrubland", "3" = "Shrubland", "4" = "Shrubland",
         "5" = "Grassland", "6" = "Grassland", "7" = "Grassland", "8" = "Grassland")

sourcetracker.meta.16S <- data.frame("SampleID" = row.names(rared.otu.16S %>% 
                                                              filter(!(grepl("LD", row.names(.)))) %>% 
                                                              rownames_to_column("SampleID") %>% as_tibble() %>%
                                                              mutate(SampleID = ifelse(SampleID == "LOT1A3_dup_new_protocol", "LOT1A3", SampleID)) %>%
                                                              column_to_rownames("SampleID"))) %>%  # just get rownames out 
  mutate(SourceSink = ifelse(grepl("Air|Soil|Env", SampleID), "source", "sink")) %>% 
  mutate(Env = str_match(SampleID,"^([a-zA-SU-Z]+)_?T?[0-9]")[,2]) %>% 
  mutate(Env = ifelse(grepl("LO", Env), 
                      key[str_match(SampleID, "[ABC]([1-8])")[,2]], # converting from place on transect to ecosystem
                      Env)) %>%
  mutate(Env = ifelse(Env == "EnvG", "Grassland_litter", Env)) %>% 
  mutate(Env = ifelse(Env == "EnvS", "Shrubland_litter", Env)) %>% 
  rename("#SampleID" = SampleID)
write.table(sourcetracker.meta.16S, file="sourcetracker.meta.16S.tsv", sep="\t")

## Repeat for fungi 
sourcetracker.OTU.ITS <- rared.otu.ITS %>% 
  .[!grepl("LD", rownames(.)), ] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("#OTU ID")
write.table(sourcetracker.OTU.ITS, file = "sourcetracker.OTU.ITS.tsv", sep="\t") ## add "# Constructed from biom file" in terminal

sourcetracker.meta.ITS <- data.frame("SampleID" = row.names(rared.otu.ITS %>% 
                                                              filter(!(grepl("LD", row.names(.)))))) %>%  
  mutate(SourceSink = ifelse(grepl("Air|Soil|Env", SampleID), "source", "sink")) %>% 
  mutate(Env = str_match(SampleID,"^([a-zA-SU-Z]+)_?T?[0-9]")[,2]) %>% 
  mutate(Env = ifelse(grepl("LO", Env), 
                      key[str_match(SampleID, "[ABC]([1-8])")[,2]], # converting from place on transect to ecosystem
                      Env)) %>%
  mutate(Env = ifelse(Env == "EnvG", "Grassland_litter", Env)) %>% 
  mutate(Env = ifelse(Env == "EnvS", "Shrubland_litter", Env)) %>% 
  rename("#SampleID" = SampleID)

write.table(sourcetracker.meta.ITS, file="sourcetracker.meta.ITS.tsv", sep="\t")

## Run scripts below via terminal
## adjust the metadata file and rared OTU table above in terminal so the column names are in the right spot for SourceTracker
## ran into an issue with the downloaded verion of sourcetracker using the installation method outlined on their github
## to fix the numpy issue, redownload sourcetracker from the developers: 
# pip install https://github.com/biota/sourcetracker2/archive/master.zip
# cd ~/Google\ Drive/My\ Drive/Professional_development/Kendra_Ch3_Manuscript/03_Analysis
# biom convert -i sourcetracker.OTU.16S.tsv -o sourcetracker.OTU.16S.biom --to-hdf5 --table-type="OTU table"
# sourcetracker2 gibbs -i sourcetracker.OTU.16S.biom -m sourcetracker.meta.16S.tsv --source_rarefaction_depth 0 --sink_rarefaction_depth 0 --jobs 2 -o sourcetracker_16S/
## 16S only took 5 minutes, do not need to rarefy again because put in rared table
# biom convert -i sourcetracker.OTU.ITS.tsv -o sourcetracker.OTU.ITS.biom --to-hdf5 --table-type="OTU table"
# sourcetracker2 gibbs -i sourcetracker.OTU.ITS.biom -m sourcetracker.meta.ITS.tsv --source_rarefaction_depth 0 --sink_rarefaction_depth 0 --jobs 2 -o sourcetracker_ITS/
## ITS took 10 minutes, do not need to rarefy again because put in rared table

## now make the sourcetracker figure
sourcetracker.proportions.16S <- as.data.frame(fread("sourcetracker_16S/mixing_proportions.txt")) %>% 
  column_to_rownames(var="V1") %>% 
  t() %>% 
  as.data.frame() %>%
  cbind(., metadata.16S.nc[rownames(.),]) %>% 
  mutate(Organism = c("Bacteria")) %>% 
  rownames_to_column(var="SampleID") %>% 
  dplyr::select(SampleID, Grassland_litter, Shrubland_litter, Soil, Unknown, Substrate, Treatment, Time_Point, Transect, Ecosystem, Location, Organism) %>% 
  pivot_longer(Grassland_litter:Unknown, names_to = "Source", values_to = "Proportion")

sourcetracker.proportions.ITS <- as.data.frame(fread("sourcetracker_ITS/mixing_proportions.txt")) %>% 
  column_to_rownames(var="V1") %>% 
  t() %>% 
  as.data.frame() %>%
  cbind(., metadata.ITS.nc[rownames(.),]) %>% 
  mutate(Organism = c("Fungi")) %>% 
  rownames_to_column(var="SampleID") %>% 
  dplyr::select(SampleID, Air, Grassland_litter, Shrubland_litter, Soil, Unknown, Substrate, Treatment, Time_Point, Transect, Ecosystem, Location, Organism) %>% 
  pivot_longer(Air:Unknown, names_to = "Source", values_to = "Proportion")

## combine
sourcetracker.proportions.combo <- rbind(sourcetracker.proportions.16S, sourcetracker.proportions.ITS)
write.table(sourcetracker.proportions.combo, file="sourcetracker_proportions_combo.tsv", sep="\t")

## summary of the average proportions explained by each source 
sourcetracker.summary <- sourcetracker.proportions.combo %>% 
  ddply(., c("Organism", "Source"), summarise,
        Avg_proportion = mean(Proportion),
        sd = sd(Proportion),
        n = n(),
        se=sd/sqrt(n))

##### Creating SourceTracker plot #####
## Create df for Dunn's mult comp test letters
library("FSA")

###### Bacteria ######
## split by ecosystem
shrubland.16S <- sourcetracker.proportions.16S %>% filter(Ecosystem == "Shrubland")
grassland.16S <- sourcetracker.proportions.16S %>% filter(Ecosystem == "Grassland")

dunnTest(Proportion ~ Source, data = shrubland.16S, method = "bonferroni")
dunnTest(Proportion ~ Source, data = grassland.16S, method = "bonferroni")

bacteria <- data.frame(Source = rep(c("Environmental_Grass", "Environmental_Shrub", "Soil", "Unknown"),2), 
                       Letters = c("A", "A", "B", "C", 
                                   "A", "A", "B", "C"), 
                       Ecosystem = c(rep("Shrubland", 4), rep("Grassland", 4)), 
                       Organism = rep("Bacteria", 8))

###### Fungi ######
shrubland.ITS <- sourcetracker.proportions.ITS %>% filter(Ecosystem == "Shrubland")
grassland.ITS <- sourcetracker.proportions.ITS %>% filter(Ecosystem == "Grassland")

dunnTest(Proportion ~ Source, data = shrubland.ITS, method = "bonferroni")
dunnTest(Proportion ~ Source, data = grassland.ITS, method = "bonferroni")

fungi <- data.frame(Source = rep(c("Air", "Environmental_Grass", "Environmental_Shrub", "Soil", "Unknown"),2), 
                    Letters = c("A", "BC", "A", "C", "AB", 
                                "A", "A", "B", "C", "C"), 
                    Ecosystem = c(rep("Shrubland", 5), rep("Grassland", 5)), 
                    Organism = rep("Fungi", 10))

## Together
letters <- bacteria %>% bind_rows(fungi)

# Make graph
ggplot(data = sourcetracker.proportions.combo) + 
  geom_text(data = letters, 
            position = position_dodge2(width = 0.8),
            aes(x = Ecosystem, y = 1.1, label = Letters)) +
  geom_boxplot(color = "black", 
               aes(x = Ecosystem, y = Proportion, fill = Source), 
               position = position_dodge2(width = 0.9)) + 
  facet_wrap(~ Organism) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_brewer(palette = "Set1", 
                    labels = set_names(c("Air", "Grass Litter", "Shrub Litter", "Soil", "Unknown"), 
                                       levels(sourcetracker.proportions.combo$Source))) +
  labs(x = "", y = "Proportion of community comp. explained", fill = "Source") +
  theme_test() + 
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=18, 
                                  margin = margin(4, 0, 4, 0)),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour="black", size = 14), 
        axis.text.y = element_text(colour="black", size = 14), 
        axis.title.x = element_text(vjust = 0.35, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 15))

ggsave(plot=last_plot(), file="sourcetracker.results.v2.pdf", width = 10, height = 5)


