## This script is to process the 16S and ITS samples from Kendra's Chapter 3 Landscape Dispersal Project
## Will use her original script "clean_glass_slide_[16S/ITS].R" as a guide for this script (showing which samples to remove vs keep)
## Output: rarefied and rounded OTU table
## Output: alpha-diversity (Shannon, richness, and Simpson) metrics for each sample in a single table
## Output: average Bray-Curtis matrix rarefied to sample depth as rarefied OTU table
## Kendra created the original script on May 20th 2021
## KB created this script on 12/12/2024 after re-processing all sequences through Qiime2 and updated it on 05/27/2025 to only include code for the manuscript

## Only the forward reads (single end) were used for downstream analysis due to poor reverse read quality (this is common)

#Reset R's Brain
rm(list=ls())

## Load packages
library(vegan)
library(plyr)
library(dplyr)
library(devtools)
library(EcolUtils)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)
library(patchwork)
library(multcomp)
library(forcats)

## set working directory
setwd("~/Google Drive/My Drive/Professional_development/Kendra_Ch3_Manuscript/03_Analysis/")

## import metadata
metadata.16S <- read.delim("KB_qiime2_redo/Qiime2_output/q1_16S_barcodes_and_pooling__metadata.tsv",  check.names = FALSE) 
metadata.ITS <- read.delim("KB_qiime2_redo/Qiime2_output/q1_ITS_barcodes_and_pooling__metadata.tsv", check.names = FALSE) %>% column_to_rownames(., var="#SampleID")

## remove controls (all negative controls and PCR mock communities look good ~ negatives are below rarefaction and mock communities match expected taxa)
## remove samples contaminated at the PCR step
metadata.16S.nc <- metadata.16S %>% filter(!Ecosystem=="NA") %>% column_to_rownames(., var="SampleID")
metadata.ITS.nc <- metadata.ITS %>% filter(!Ecosystem=="NA")

## import 16S and ITS OTU table
otu.table.16S <- read.delim("KB_qiime2_redo/Qiime2_output/16S_table_filtered.tsv", check.names = FALSE) %>% column_to_rownames(., var = "OTU_ID")
otu.table.ITS <- read.delim("KB_qiime2_redo/Qiime2_output/ITS_table_filtered.tsv", check.names = FALSE) %>% column_to_rownames(., var = "OTU_ID")

## Remove PCR + and - (checked all negatives are below rarefaction depth), also remove samples that Kendra originally removed due to errors 
otu.table.16S <- otu.table.16S %>% 
  .[, !grepl("PCR|Neg|AirT4", colnames(.))] %>% # AirT4 samples fell on the ground and are not truly air samples
  dplyr::select(-c(LDT2R10CG, LOT4A1, LDT0R5G_oops)) # these two samples were contaminated in the PCR step, LDT0R5G_oops accidentally added 10 uL

otu.table.ITS <- otu.table.ITS %>% 
  .[, !grepl("PCR|Neg|AirT4", colnames(.))] # AirT4 samples fell on the ground and are not truly air samples

## import taxonomy 
## everything assigned to with confidence >70% at the speacies level (will never use the species level assignment)
taxa.table.16S <- read.delim("KB_qiime2_redo/Qiime2_output/16S_taxonomy.tsv", check.names = FALSE)
taxa.table.ITS <- read.delim("KB_qiime2_redo/Qiime2_output/ITS_taxonomy.tsv", check.names = FALSE)

## get taxa levels for each OTU, parse by taxonomic rank
## divide column using ";" and convert list to dataframe
taxa.table.16S.list  <- ldply(str_split(taxa.table.16S$Taxon, pattern = ";"), rbind)
names(taxa.table.16S.list) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa.table.16S.list <- lapply(taxa.table.16S.list, function(x){gsub(pattern = "[a-z]__", "", x)}) %>% 
  as.data.frame()
## get final table with subdivided information
taxa.table.16S.final <- cbind(taxa.table.16S[,1:2],taxa.table.16S.list) %>% 
  column_to_rownames(., var = 'Feature ID')

## now do the same for ITS
## get taxa levels for each OTU, parse by taxonomic rank
## divide column using ";" and convert list to dataframe
taxa.table.ITS.list  <- ldply(str_split(taxa.table.ITS$Taxon, pattern = ";"), rbind)
names(taxa.table.ITS.list) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_hypothesis")
taxa.table.ITS.list <- lapply(taxa.table.ITS.list, function(x){gsub(pattern = "[a-z]__", "", x)}) %>% as.data.frame()
## remove the remaining s from the species hypothesis row
taxa.table.ITS.list$Species_hypothesis <- gsub("s", "", taxa.table.ITS.list$Species_hypothesis)
## get final table with subdivided information
taxa.table.ITS.final <- cbind(taxa.table.ITS[,1:2],taxa.table.ITS.list) %>% 
  column_to_rownames(., var = 'Feature ID')

## Merge taxonomy to OTU IDs, filter out desired domain if needed, remove unnecessary columns, and change rownames to OTU ID
## For the 16S table, there are some OTUs associated with Archea or Eukaryotes so remove those (rowsums = 0)
otu.tax.table.16S <- as.data.frame(merge(taxa.table.16S.final, otu.table.16S, by.x = "row.names", by.y = "row.names")) %>% 
  filter(Domain=="Bacteria") %>% 
  column_to_rownames(., var = "Row.names") %>% 
  .[,-1] %>%  
  mutate(rowsums = rowSums(.[,-c(1:8)])) %>% 
  filter(rowsums!=0) %>% subset(., select = -rowsums) %>% 
  as.data.frame()

## save the table
write.table(otu.tax.table.16S, file = "otu.tax.table.16S.v1.tsv", quote = FALSE, sep = "\t")
## import the table if coming back in
#otu.tax.table.16S <- read.delim("otu.tax.table.16S.v1.tsv", sep="\t", check.names = FALSE)

otu.tax.table.ITS <- as.data.frame(merge(taxa.table.ITS.final, otu.table.ITS, by.x = "row.names", by.y = "row.names")) %>% 
  column_to_rownames(., var = "Row.names") %>% 
  .[,-1]

## save the table 
write.table(otu.tax.table.ITS, file = "otu.tax.table.ITS.v1.tsv", quote = FALSE, sep = "\t")
## import the table if coming back in
#otu.tax.table.ITS <- read.delim("otu.tax.table.ITS.v1.tsv", sep="\t", check.names = FALSE)

## now remove the taxonomic information to look at read distribution across samples
otu.table.clean.16S <- otu.tax.table.16S[,-c(1:7)] %>% as.data.frame() #17755 OTUs (or ESVs), 171 samples
write.table(otu.table.clean.16S, file = "otu.table.clean.16S.v1.tsv", quote = FALSE, sep = "\t")
# otu.table.clean.16S <- read.delim("otu.table.clean.16S.v1.tsv", sep="\t")

otu.table.clean.ITS <- otu.tax.table.ITS[,-c(1:8)] %>% as.data.frame() #8408 OTUs (or ESVs), 162 samples
write.table(otu.table.clean.ITS, file = "otu.table.clean.ITS.v1.tsv", quote = FALSE, sep = "\t")
# otu.table.clean.ITS <- read.delim("otu.table.clean.ITS.v1.tsv", sep="\t")

##### Rarefaction depth #####
## Rarefying - normalizes read depth across all samples. 
## Allows for an equal comparison across different sample types, at the risk of excluding rarer taxa
## A "good" rarefaction depth should minimize sample loss while maximizing OTU richness.

## first look at read depth across samples
###### 16S ######
hist(colSums(otu.table.clean.16S))
sort(colSums(otu.table.clean.16S)) ## all negatives below 970, highest negative = 876 (highest ITS negative = 547)
## a lot of low sequence #s because low abundance samples
sums <- colSums(otu.table.clean.16S)
length(sums[sums >=970]) / length(sums) ## rarefying to 970 would keep 66% of samples (Kendra originally rarefied to 971 which kept 63.8% of samples)
mean(sums) #16992
min(sums) #0 (all eukaryotic)
max(sums) #116791
median(sums) #4483

## Will rarefy at 970
barplot(sort(colSums(otu.table.clean.16S)), ylim=c(0, max(colSums(otu.table.clean.16S))),
        xlim = c(0, NCOL(otu.table.clean.16S)), col = "Gray", ylab = "Read Depth", xlab = "Sample")
unrarefied.cruves <- rarecurve(t(otu.table.clean.16S), step = 1000, label = F)
rared.curves <- rarecurve(t(otu.table.clean.16S), step = 1000, label = F, sample = 970)
## will rarefy 16S table to 970 (will lose 57 samples)

###### ITS ######
hist(colSums(otu.table.clean.ITS))
sort(colSums(otu.table.clean.ITS)) ## all negatives below 1840 , highest ITS negative = 547
## a lot of low sequence #s because low abundance samples
sums.ITS <- colSums(otu.table.clean.ITS)
length(sums.ITS[sums.ITS >=1840 ]) / length(sums.ITS) ## rarefying to 1840 would keep 74% of samples (Kendra originally rarefied to 2400 which kept 80% of samples)
mean(sums.ITS) #35758
min(sums.ITS) #1 *all eukaryotic
max(sums.ITS) #175960
median(sums.ITS) #30010

barplot(sort(colSums(otu.table.clean.ITS)), ylim=c(0, max(colSums(otu.table.clean.ITS))),
        xlim = c(0, NCOL(otu.table.clean.ITS)), col = "Gray", ylab = "Read Depth", xlab = "Sample")
unrarefied.cruves <- vegan::rarecurve(t(otu.table.clean.ITS), step = 1000, label = FALSE)
rared.curves <- rarecurve(t(otu.table.clean.ITS), step = 1000, sample = 8000, label = FALSE, col = "blue")
## will rarefy ITS table to 1840 (will lose 42 samples)

##### Beta diversity (Fig 2B) #####
## beta diversity is an assessment of community diversity between samples
## Bray-Curtis dissimilarity is an index that measures how different (or similar) communities from two samples are from each other
## Bray-Curtis is calculated as (2*(# OTUs shared between the Sample A and B))/(total OTU in Sample A and total OTU in Sample B)
set.seed(100)
###### 16S ######
otu.table.clean.16S.t <- t(otu.table.clean.16S) %>% as.data.frame()
bray.rared.16S <- as.matrix(avgdist(otu.table.clean.16S.t, 970, transf = sqrt, iterations = 1000,
                                    meanfun = median, dmethod = "bray"))
## double check that duplicated samples only have one representative, all samples are unique
write.table(bray.rared.16S, file = "bray.rared.16S.tsv", quote = FALSE, sep = "\t")
# bray.rared.16S <- read.delim("bray.rared.16S.tsv")

## save table for PRIMER analysis (where PERMANOVA statistics are conducted)
metadata.16S.nc.t <- as.data.frame(t(metadata.16S.nc))
bray.rared.16S.meta <- rbind(bray.rared.16S, metadata.16S.nc.t[, colnames(bray.rared.16S)])
write.table(bray.rared.16S.meta, "bray.rared.16S.meta.tsv", sep="\t")

## preform NMDS (non-metric multidimensional scaling)
NMDS.16S <- metaMDS(bray.rared.16S, distance="bray", k=2, engine="monoMDS", maxit=500, try=80, trymax=100)
NMDS.16S.points <- as.data.frame(NMDS.16S[["points"]]) %>% cbind(., metadata.16S.nc[rownames(.),])
write.table(NMDS.16S.points, file="NMDS.16S.points.v1.tsv", sep="\t")
# NMDS.16S.points <- read.delim("NMDS.16S.points.v1.tsv")

###### ITS ######
otu.table.clean.ITS.t <- t(otu.table.clean.ITS) %>% as.data.frame()
bray.rared.ITS <- as.matrix(avgdist(otu.table.clean.ITS.t, 1840, transf = sqrt, iterations = 1000,
                                    meanfun = median, dmethod = "bray")) ##no sample errors
write.table(bray.rared.ITS, file = "bray.rared.ITS.tsv", quote = FALSE, sep = "\t")
# bray.rared.ITS <- read.delim("bray.rared.ITS.tsv")

## save table for PRIMER analysis (where PERMANOVA statistics are conducted)
metadata.ITS.nc.t <- as.data.frame(t(metadata.ITS.nc))
bray.rared.ITS.meta <- rbind(bray.rared.ITS, metadata.ITS.nc.t[, colnames(bray.rared.ITS)])
write.table(bray.rared.ITS.meta, "bray.rared.ITS.meta.tsv", sep="\t")

## preform NMDS (non-metric multidimensional scaling)
NMDS.ITS <- metaMDS(bray.rared.ITS, distance="bray", k=2, engine="monoMDS", maxit=500, try=80, trymax=100) ## high stress, over 0.25
NMDS.ITS.points <- as.data.frame(NMDS.ITS[["points"]]) %>% cbind(., metadata.ITS[rownames(.),])
write.table(NMDS.ITS.points, "NMDS.ITS.points.v1.tsv", sep = "\t")
# NMDS.ITS.points <- read.delim("NMDS.ITS.points.v1.tsv")

###### NMDS Plots #####
nmds.theme <- theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1, linetype = "solid"),
                    axis.line = element_blank(),
                    panel.background = element_blank(),
                    axis.text = element_text(size = 12, colour = "black"),
                    axis.title = element_text(size = 12),
                    legend.title = element_text(size = 12),
                    legend.text = element_text(size = 11),
                    legend.background = element_blank(),
                    legend.box.background = element_blank(),
                    legend.key = element_blank(),
                    axis.ticks.length = unit(.25, "cm"))

### subset the open glass slides only
bray.rared.16S.open <- bray.rared.16S %>% as.data.frame() %>% 
  filter(grepl("LO", rownames(.))) %>% 
  .[, grep("LO", names(.))]

bray.rared.16S.open.meta <- rbind(bray.rared.16S.open, metadata.16S.nc.t[, colnames(bray.rared.16S.open)])
write.table(bray.rared.16S.open.meta, "bray.rared.16S.open.meta.tsv", sep="\t") ## save for PRIMER

### subset the open glass slides only
bray.rared.ITS.open <- bray.rared.ITS %>% as.data.frame() %>% 
  filter(grepl("LO", rownames(.))) %>% 
  .[, grep("LO", names(.))]

bray.rared.ITS.open.meta <- rbind(bray.rared.ITS.open, metadata.ITS.nc.t[, colnames(bray.rared.ITS.open)])
write.table(bray.rared.ITS.open.meta, "bray.rared.ITS.open.meta.tsv", sep="\t") ## save for PRIMER

## subset death slides
bray.rared.16S.death <- bray.rared.16S %>% as.data.frame() %>% 
  filter(grepl("LD", rownames(.))) %>% 
  .[, grep("LD", names(.))]

bray.rared.16S.death.meta <- rbind(bray.rared.16S.death, metadata.16S.nc.t[, colnames(bray.rared.16S.death)])
write.table(bray.rared.16S.death.meta, "bray.rared.16S.death.meta.tsv", sep="\t") ## save for PRIMER

## subset death slides
bray.rared.ITS.death <- bray.rared.ITS %>% as.data.frame() %>% 
  filter(grepl("LD", rownames(.))) %>% 
  .[, grep("LD", names(.))]

bray.rared.ITS.death.meta <- rbind(bray.rared.ITS.death, metadata.ITS.nc.t[, colnames(bray.rared.ITS.death)])
write.table(bray.rared.ITS.death.meta, "bray.rared.ITS.death.meta.tsv", sep="\t") ## save for PRIMER

## NMDS for accumulation (open) slides only
## 16S
NMDS.16S.open <- metaMDS(bray.rared.16S.open, distance="bray", k=2, engine="monoMDS", maxit=500, try=80, trymax=100)
NMDS.16S.open.points <- as.data.frame(NMDS.16S.open[["points"]]) %>% 
  cbind(., metadata.16S.nc[rownames(.),]) %>% 
  mutate(Organism = c("Bacteria")) %>% 
  rownames_to_column(var="SampleID") %>% 
  dplyr::select(SampleID, MDS1, MDS2, Substrate, Treatment, Time_Point, Transect, Ecosystem, Location, Organism, Distance_to_midline)

write.table(NMDS.16S.open.points, file="NMDS.16S.open.points.v1.tsv", sep="\t")
# NMDS.16S.open.points <- read.delim("NMDS.16S.open.points.v1.tsv")

## ITS
NMDS.ITS.open <- metaMDS(bray.rared.ITS.open, distance="bray", k=2, engine="monoMDS", maxit=500, try=80, trymax=100)
NMDS.ITS.open.points <- as.data.frame(NMDS.ITS.open[["points"]]) %>% 
  cbind(., metadata.ITS.nc[rownames(.),]) %>% 
  mutate(Organism = c("Fungi")) %>% 
  rownames_to_column(var="SampleID") %>% 
  dplyr::select(SampleID, MDS1, MDS2, Substrate, Treatment, Time_Point, Transect, Ecosystem, Location, Organism, Distance_to_midline)

write.table(NMDS.ITS.open.points, file="NMDS.ITS.open.points.v1.tsv", sep="\t")
# NMDS.ITS.open.points <- read.delim("NMDS.ITS.open.points.v1.tsv")

NMDS.meta.combined <- rbind(NMDS.16S.open.points, NMDS.ITS.open.points)

## Make Fig 2B
##below is copied from Kendra's original scripts
# Create our visualization of dispersal slides + sources
cols <- c("#508104", "#db8200") # slightly darker set of orange + green

NMDS.combo.ecosystem.plot <- ggplot(data = NMDS.meta.combined) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Ecosystem)), size = 4, alpha = 0.8, stroke=1) +
  facet_wrap(~ Organism, 
             scales = "free") + 
  theme_test() +
  scale_color_manual(values = cols, name = "Ecosystem") +
  labs(x="NMDS1", y="NMDS2") +
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=17, face = "bold", 
                                  margin = margin(4, 0, 4, 0)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(vjust = 0.35, size = 12),
        axis.title.y = element_text(vjust = 0.35, size = 12), 
        legend.title = element_text(face = "bold", size = 15), 
        legend.text = element_text(size = 12)) 

plot(NMDS.combo.ecosystem.plot)
ggsave(plot=last_plot(), file="NMDS.dispersal.communities.by.ecosystem.bact.fungi.pdf")
#pdf("../../04_Figures/NMDS_Dispersal_Communities_by_Ecosystem_Fungi_Bacteria.pdf", width = 8, height = 4)


##### Mantel Test (Fig 2C) #####
# will use Kendra's pre-processed matrices for plant community composition around each sample and distance between samples
# will need to remake the bray curtis dissimilarity matrix based on the new OTU table from qiime2
library(ncf)

###### Bacteria #####

## Microbial matrix
M.microbial.16S <- bray.rared.16S.open %>% 
  rename(LOT1A3 = LOT1A3_dup_new_protocol) %>% 
  rownames_to_column("SampleID") %>% 
  mutate(SampleID = ifelse(SampleID == "LOT1A3_dup_new_protocol", "LOT1A3", SampleID)) %>% 
  column_to_rownames("SampleID") %>% as.matrix

## Distance Matrix
# Load site coordinates and use to make dataframe of coordinates for each sample
site.coords <- as.data.frame(fread("Kendra_processed_data/site_locs_coordinates.tsv")) %>%
  rename(Site = V1)

M.distance.16S <- data.frame("Sample" = row.names(M.microbial.16S)) %>% # start with same names in same order
  mutate(Site = str_extract(Sample, "[ABC][1-8]")) %>%  # get site ID
  left_join(site.coords) %>% # join with our coordinates
  column_to_rownames("Sample") %>%
  dplyr::select(x,y) %>% dist(., method = "euclidean") %>% # create Euclidean distance matrix
  as.matrix

# mantel time!! this is testing correlation between microbial composition + distance
distance.mantel.output <- vegan::mantel(M.microbial.16S, M.distance.16S, permutations = 10000, method = "pearson")
distance.mantel.output$statistic # r = 0.06276195
distance.mantel.output$signif # p = 0.1269873, bacterial composition is not significantly correlated with distance between samples

# 1 meter PLANT MATRIX
## From Kendra's original script "For this analysis, we are arbitrarily using the plant composition within 1 meter of each site. We could have chosen a different radius here, but this one seemed OK."
plants.1m <- fread("Kendra_processed_data/map/grass_shrub_1_meter_composition_df.tsv", data.table = FALSE) %>% 
  rename(Site = V1)

M.plants.1m.16S <- data.frame("Sample" = row.names(M.microbial.16S)) %>% 
  mutate(Site = str_extract(Sample, "[ABC][1-8]")) %>% 
  left_join(plants.1m) %>% 
  column_to_rownames("Sample") %>% 
  dplyr::select(!(Site)) %>% 
  vegdist() %>% as.matrix

# mantel time!! this is testing correlation between microbial composition + plant composition @ 1m 
plant.1m.mantel.output <- vegan::mantel(M.microbial.16S, M.plants.1m.16S, permutations = 10000, method = "pearson")
plant.1m.mantel.output$statistic # r = 0.1245136
plant.1m.mantel.output$signif # p = 0.04129587

## partial mantel ~ redo mantel tests while controlling for a third matrix (don't need to do because bacterial comp is not sig. impacted by distance between samples)
partial.mantel.plant1.output <- partial.mantel.test(M.microbial.16S, M.plants.1m.16S, M.distance.16S, resamp = 1000, method = "pearson", quiet = TRUE)
partial.mantel.plant1.output$MantelR
partial.mantel.plant1.output$p

### Here we start the loop through our difference distances

## PLANT MATRIX
## reading in plant composition bc and make df of plant composition of each sample
plant.list <- list.files("Kendra_processed_data/map/", pattern = "_meter_composition_df.tsv")
r.list <- list()
p.list <- list()
radius.list <- list()
partial.r.list <- list()
partial.p.list <- list()

## The radius at which you measured plant composition
radius <- c(0.1, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4) # COPY THIS FROM THE -- MAP_make_plant_composition_spatial_polygon_file.R

for (i in 1:length(plant.list)) {
  plants <- fread(paste0("Kendra_processed_data/map/", plant.list[[i]]), data.table = FALSE) %>% 
    rename(Site = V1)
  
  M.plants.16S <- data.frame("Sample" = row.names(M.microbial.16S)) %>% 
    mutate(Site = str_extract(Sample, "[ABC][1-8]")) %>% 
    left_join(plants) %>% 
    column_to_rownames("Sample") %>% 
    dplyr::select(!(Site)) %>% 
    vegdist() %>% as.matrix
  
  # mantel time!! this is testing correlation between microbial composition + plant composition
  # mantel.output <- vegan::mantel(M.microbial.16S, M.plants.16S, permutations = 10000, method = "pearson")
  # 
  # 
  # r.list[[i]] <- mantel.output$statistic
  # p.list[[i]] <- mantel.output$signif
  
  ## Now do a partial Mantel, select for test for microbe x plant (where variation from distance matrix is removed (r13)
  partial.mantel.output <- partial.mantel.test(M.microbial.16S, M.distance.16S, M.plants.16S, resamp = 10000, method = "pearson", quiet = TRUE)
  # print(radius[i])
  # print(partial.mantel.output)
  #
  partial.r.list[[i]] <- partial.mantel.output$MantelR[["r13"]]
  partial.p.list[[i]] <- partial.mantel.output$p[[2]]
  
}
# mantel.corr.all.16S <- data.frame("r" = unlist(r.list), 
#                                             "p" = unlist(p.list), 
#                                             "radius" = radius) %>% 
#   mutate(signif = ifelse(p < 0.05, "yes", "no"), 
#          Samples = "all", 
#          Organism = c("Bacteria")) ##radius 1.25 is not correlated anymore, p = 0.052 so not being included as significant (was for Kendra)

partial.mantel.corr.all.16S <- data.frame("r" = unlist(partial.r.list), 
                                          "p" = unlist(partial.p.list), 
                                          "radius" = radius) %>% 
  mutate(signif = ifelse(p < 0.05, "yes", "no"), 
         Samples = "all", 
         Organism = c("Bacteria"))

###### Fungi #####

## Microbial matrix
M.microbial.ITS <- bray.rared.ITS.open %>% 
  as.matrix

## Distance Matrix
# Load site coordinates and use to make dataframe of coordinates for each sample
site.coords <- as.data.frame(fread("Kendra_processed_data/site_locs_coordinates.tsv")) %>%
  rename(Site = V1)

M.distance.ITS <- data.frame("Sample" = row.names(M.microbial.ITS)) %>% # start with same names in same order
  mutate(Site = str_extract(Sample, "[ABC][1-8]")) %>%  # get site ID
  left_join(site.coords) %>% # join with our coordinates
  column_to_rownames("Sample") %>%
  dplyr::select(x,y) %>% dist(., method = "euclidean") %>% # create Euclidean distance matrix
  as.matrix

# mantel time!! this is testing correlation between microbial composition + distance
distance.mantel.output.ITS <- vegan::mantel(M.microbial.ITS, M.distance.ITS, permutations = 10000, method = "pearson")
distance.mantel.output.ITS$statistic # r = 0.1898438
distance.mantel.output.ITS$signif # p = 0.00029997 ## fungal composition is significantly impacted by distance between samples

# 1 meter PLANT MATRIX
## From Kendra's original script "For this analysis, we are arbitrarily using the plant composition within 1 meter of each site. We could have chosen a different radius here, but this one seemed OK."

plants.1m <- fread("Kendra_processed_data/map/grass_shrub_1_meter_composition_df.tsv", data.table = FALSE) %>% 
  rename(Site = V1)

M.plants.1m.ITS <- data.frame("Sample" = row.names(M.microbial.ITS)) %>% 
  mutate(Site = str_extract(Sample, "[ABC][1-8]")) %>% 
  left_join(plants.1m) %>% 
  column_to_rownames("Sample") %>% 
  dplyr::select(!(Site)) %>% 
  vegdist() %>% as.matrix

# mantel time!! this is testing correlation between microbial composition + plant composition @ 1m 
plant.1m.mantel.output.ITS <- vegan::mantel(M.microbial.ITS, M.plants.1m.ITS, permutations = 10000, method = "pearson")
plant.1m.mantel.output.ITS$statistic # r = 0.6073828
plant.1m.mantel.output.ITS$signif # p = 9.999e-05

## partial mantel ~ redo mantel tests while controlling for a third matrix
partial.mantel.plant1.output.ITS <- partial.mantel.test(M.microbial.ITS, M.distance.ITS, M.plants.1m.ITS, resamp = 1000, method = "pearson", quiet = TRUE)

### Here we start the loop through our difference distances

## PLANT MATRIX
## reading in plant composition bc and make df of plant composition of each sample
plant.list <- list.files("Kendra_processed_data/map/", pattern = "_meter_composition_df.tsv")
r.list.ITS <- list()
p.list.ITS <- list()
radius.list <- list()
partial.r.list.ITS <- list()
partial.p.list.ITS <- list()

## The radius at which you measured plant composition
radius <- c(0.1, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4) # COPY THIS FROM THE -- MAP_make_plant_composition_spatial_polygon_file.R

for (i in 1:length(plant.list)) {
  plants <- fread(paste0("Kendra_processed_data/map/", plant.list[[i]]), data.table = FALSE) %>% 
    rename(Site = V1)
  
  M.plants.ITS <- data.frame("Sample" = row.names(M.microbial.ITS)) %>% 
    mutate(Site = str_extract(Sample, "[ABC][1-8]")) %>% 
    left_join(plants) %>% 
    column_to_rownames("Sample") %>% 
    dplyr::select(!(Site)) %>% 
    vegdist() %>% as.matrix
  
  # # mantel time!! this is testing correlation between microbial composition + plant composition
  # mantel.output.ITS <- vegan::mantel(M.microbial.ITS, M.plants.ITS, permutations = 10000, method = "pearson")
  # 
  # 
  # r.list.ITS[[i]] <- mantel.output.ITS$statistic
  # p.list.ITS[[i]] <- mantel.output.ITS$signif
  
  ## Now do a partial Mantel, select for test for microbe x plant where variation from distance matrix is removed (r13.2; sig. correlated with distance matrix)
  partial.mantel.output.ITS <- partial.mantel.test(M.microbial.ITS, M.distance.ITS, M.plants.ITS, resamp = 10000, method = "pearson", quiet = TRUE)
  print(radius[i])
  print(partial.mantel.output.ITS)
  
  partial.r.list.ITS[[i]] <- partial.mantel.output.ITS$MantelR[["r13.2"]]
  partial.p.list.ITS[[i]] <- partial.mantel.output.ITS$p[[5]]
  
}

# mantel.corr.all.ITS <- data.frame("r" = unlist(r.list.ITS), 
#                                   "p" = unlist(p.list.ITS), 
#                                   "radius" = radius) %>% 
#   mutate(signif = ifelse(p < 0.05, "yes", "no"), 
#          Samples = "all", 
#          Organism = c("Fungi")) ##all significantly correlated

partial.mantel.corr.all.ITS <- data.frame("r" = unlist(partial.r.list.ITS), 
                                          "p" = unlist(partial.p.list.ITS), 
                                          "radius" = radius) %>% 
  mutate(signif = ifelse(p < 0.05, "yes", "no"), 
         Samples = "all", 
         Organism = c("Fungi"))

## combine
# mantel.corr.all.combo <- mantel.corr.all.16S %>% bind_rows(mantel.corr.all.ITS)
partial.mantel.corr.all.combo <- partial.mantel.corr.all.16S %>% bind_rows(partial.mantel.corr.all.ITS)

## plotting copied from Kendra's original script

## GRAPH TIME
mantel.combo.plot <- ggplot(data = partial.mantel.corr.all.combo, aes(x = radius, y = r)) +
  geom_line(aes(linetype = Samples), show.legend = FALSE) + 
  geom_point(color = "white", size = 5) +
  geom_point(aes(shape = signif), 
             colour = "black", size = 5) +
  facet_wrap(~ Organism, scales = "free_y") +
  theme_test() + 
  xlab("Radius (m)") +
  ylab("Correlation Coefficient (r)") + 
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=17, face = "bold", 
                                  margin = margin(4, 0, 4, 0)),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour="black", size = 12), 
        axis.text.y = element_text(colour="black", size = 12), 
        axis.title.x = element_text(vjust = 0.35, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15)) +
  scale_shape_manual(values=c(1,19), labels = c("P > 0.05", "P < 0.05")) # , guide = "none"

plot(mantel.combo.plot)
ggsave(plot=last_plot(), file="Partial.mantel.correlations.at.increasing.radii.pdf")
#pdf("../../04_Figures/Mantel_Correlations_at_Increasing_Radii.pdf", width = 8, height = 4)

## make Figure 2b and c
NMDS.combo.ecosystem.plot + mantel.combo.plot + plot_layout(ncol = 1)
ggsave(plot=last_plot(), file="NMDS.partial.mantel.correlations.combo.pdf")

##### Immigration and Death rates (Fig 2A) #####
## All glass slide abundance analyses were conducted in the "clean flow cytometry data KB" and "Fig immigration and death rates" scripts
