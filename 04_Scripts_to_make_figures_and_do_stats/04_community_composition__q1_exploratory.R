# Statistics used in bacterial community analysis for leaf litter dispersal bags (Q2)
# Kendra Walters
# November 10th 2020



# Room of ()ments
require(vegan)
require(ggplot2)
require(data.table)
require(anchors)
require(EcolUtils)
require(car)
# library(biomformat)


# Setting up and reading in
setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_Data/06_Community_Composition/06_Output_from_qiime/Q1_single/")

#### Background_rarefying__dont_run -----------------------------------------------------
q1.table <- as.data.frame(fread("12_table_filtered.tsv")) 
#   --p-trim-left 5 \
#   --p-trunc-len 224 \
q1.taxonomy <- as.data.frame(fread("13_taxonomy_filtered.tsv"))
q1.combined <- as.data.frame(fread("14_table_filtered_with_taxonomy_filtered.tsv"))


sums <- colSums(q1.table[, 2:ncol(q1.table)])
hist(sums[sums < 20000], breaks = 100)

sort(sums) # maybe rarefy to 971?


# Make the OTU ID the row name
row.names(q1.table) <- q1.table$`#OTU ID`
row.names(q1.combined) <- q1.combined$`#OTU ID`

q1.table <- q1.table[ , !(names(q1.table) %in% c("#OTU ID"))]
q1.combined <- q1.combined[ , !(names(q1.combined) %in% c("#OTU ID"))]

# and transform the dataframes because vegan expects rows = samples and columns = species
q1.table <- as.data.frame(t(q1.table))
q1.combined <- as.data.frame(t(q1.combined))

# Rarefaction depth
X <- 971

drop <- row.names(q1.table[rowSums(q1.table) < X, ])

# This way will round out the averages, so you won't have exactly the same # of sequences in each, but you also won't be keeping your super rare, super low abundance << 1 taxa
q1.table.rarefied <- as.data.frame((rrarefy.perm(q1.table, sample = X, n = 1000, round.out = TRUE)))

q1.table.rarefied <- q1.table.rarefied[!(row.names(q1.table.rarefied) %in% drop), ] # drop those samples that had under rarefaction threshold
sort(rowSums(q1.table.rarefied))


write.table(q1.table.rarefied, file = "Q1_16S_single_filtered_rarefied__ROUNDED.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)



#First, filter out mock and negative community standards by name.
mock_and_neg <- c("Negative2", "PCR Mock1", "PCR Mock2", "PCR Neg1", "PCR Neg2", "Negative1")
q1.table <- q1.table[,!(names(q1.table) %in% mock_and_neg)]
q1.combined <- q1.combined[,!(names(q1.combined) %in% mock_and_neg)]








#####

##### Loading rarefied table + combining it with metadata ####

q1.table.rarefied <- as.data.frame(fread("Q1_16S_single_filtered_rarefied__ROUNDED.tsv"))
names(q1.table.rarefied) <- c("#OTUID", names(q1.table.rarefied[2:ncol(q1.table.rarefied)]))

metadata <- as.data.frame(fread("./../../04_final_pre_sequence_information/q1_16S_barcodes_and_pooling__metadata.tsv"))
metadata <- metadata[ , c("SampleID", "Substrate", "Treatment", "Time_Point", "Transect",
                          "Ecosystem", "Location")]

q1.table.rarefied.meta <- merge(q1.table.rarefied, metadata, by.x = "#OTUID", by.y = "SampleID", all.x = TRUE, all.y = FALSE)










###############---------- Alpha diversity---------  ##############
#Alpha diversity is a measure of species richness within an environment.

q1.table.rarefied <- as.data.frame(fread("OTU_filtered_plus_rarefied__ROUNDED.tsv"))
row.names(q1.table.rarefied) <- q1.table.rarefied$V1
q1.table.rarefied$V1 <- NULL

## what percent of the samples did not sequence
pooling <- as.data.frame(fread("~/Downloads/All Samples with Pooling Information - Sheet1.tsv"))
pooling <- pooling[grepl("^Air|^L", pooling$Input_ID), ]

yes.sequenced <- row.names(q1.table.rarefied)

pooling$Sequenced <- ifelse(pooling$Input_ID %in% yes.sequenced, "yes", "no")

table(pooling[pooling$`Pool Volume (uL)` == "8", ]$Sequenced)
table(pooling$`Pool Volume (uL)`)

NO_pro_K_list <- c("LUT4R3", "LOT4R5", "LUT4R5", "LOT5R6", "LUT3R4", "LUT4R4", "LET4R4", "Negative1", 
                   "LDT1R1", "LUT5R4", "LUT1R1", "LOT4R2", "LOT4R6", "LOT5R1", "LUT5R5", "LDT3R7", "LOT4R7", 
                   "LDT4R7", "Negative2", "LDT4R1", "LUT5R3", "LDT4R5", "LOT5R2", "LOT4R3", "LOT1R1", "LOT5R5", 
                   "LOT4R4", "LUT4R2", "LOT1R3", "LDT3R6", "LUT1R4", "LDT4R3", "LUT5R6", "LUT3R2")
pooling$Proteinase_K <- ifelse(pooling$Input_ID %in% NO_pro_K_list, "No",
                               ifelse(grepl("Litter", pooling$Input_ID), "No","Yes"))

ggplot(data = pooling[!(grepl("Litter", pooling$Input_ID)), ]) + geom_boxplot(aes(x = Proteinase_K, y = `Pool Volume (uL)`))

# Read in metadata!
metadata <- as.data.frame(fread("/Users/walters_kendra/Google Drive/Dispersal/Data/Community Composition Analysis/07_Bioinformatics/metadata_q1.txt"))

# Calculate alpha diversity!
shannon <- as.data.frame(diversity(q1.table.rarefied, index = "shannon"))
names(shannon) <- c("Shannon")
richness <- as.data.frame(apply(q1.table.rarefied[,]>0, 1, sum))
names(richness) <- c("Richness")
simpson <- as.data.frame(diversity(q1.table.rarefied, index = "simpson"))
names(simpson) <- c("Simpson")

# And merge them all together!
tmp1 <- merge(richness, shannon, by.x = "row.names", by.y = "row.names")
all.alpha <- merge(tmp1, simpson, by.x = "Row.names", by.y = "row.names")
names(all.alpha) <- c("SampleID", names(all.alpha)[2:length(names(all.alpha))])

#Merge with metadata to create a plot.
merged_alpha <- merge(all.alpha, metadata, by.x = "SampleID", by.y = "SampleID")




#### BETA DIVERSITY ########

###______ Making NMDS:



# First, calculate the BC matrix so you can do your rarefactions and sqrt transformations
bray.dist <- avgdist(q1.table, sample = 971, meanfun = median, transf = sqrt, iterations = 999) #default is Bray-Curtis calculated by vegdist()
bray.dist.df <- as.data.frame(as.matrix(bray.dist))
write.table(x = bray.dist.df, "BC_median__r971_single.tsv", sep = "\t", row.names = FALSE)

#Next, we will run an NMDS, which is a form of ordination and can be used to visualize beta diversity.
NMDS1 <- metaMDS(bray.dist, autotransform = FALSE, k = 2)

#This will make the first two columns as the x and y coordinates, so that they may be plotted.
coordinates <- data.frame(NMDS1$points[,1:2])

#To make a more sophisticated plot, we will merge the stress scores with the metadata.
nmds_plus_metadata <- merge(coordinates, metadata, by.x = "row.names", by.y = "SampleID")
write.table(x = nmds_plus_metadata, file = "metadata_for_PRIMER_because_stupid_sample_order_q1_single_r971.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


#Time to make a new plot colored with metadata.
nmds_plus_metadata <- as.data.frame(fread("/Users/walters_kendra/Google Drive/Dispersal/Data/Community Composition Analysis/08_Statistics/Q1_16S/second_sequencing_run/single/"))

cols <- c("Elevated (Air)"="red3",  "Timepoint 2"="darkorange2",  "Timepoint 3"= "gold", "Timepoint 4" = "green4", "Timepoint 5" = "mediumblue")

nmds_plus_metadata$Dispersal_Route <- factor(nmds_plus_metadata$Dispersal_Route, levels = c(
  "Elevated (Air)", "Overhead (Air + Vegetation)", "Open (Air + Vegetation + Soil)", 
  "Environmental Air", "Environmental Litter", "Environmental Soil"))

nmds_plus_metadata$Ecosystem_Location <- paste0(nmds_plus_metadata$Ecosystem, "_", nmds_plus_metadata$Location)


not.death <- nmds_plus_metadata[nmds_plus_metadata$Treatment != "Death", ]
not.death.and.not.env <- not.death[!(grepl("Environmental", not.death$Treatment)), ]

ggplot(data = not.death.and.not.env) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Transect)), size = 4) + #Creates and colors legend to match, modify after $ here.
  labs(col =  "Dispersal Route") +
  theme_classic() +
  scale_color_brewer(palette = "Set1") 


#Time to make a new plot colored with metadata.
cols <- c("Timepoint 1"="red3",  "Timepoint 2"="darkorange2",  "Timepoint 3"= "gold", "Timepoint 4" = "green4", "Timepoint 5" = "mediumblue")


ggplot(data = nmds_plus_metadata) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(nmds_plus_metadata$Timepoint), shape = as.factor(nmds_plus_metadata$Dispersal_Route)), size = 4) + #Creates and colors legend to match, modify after $ here.
  labs(col =  "Timepoint", shape = "Dispersal Route") +
  theme_classic() + 
  scale_color_manual(values = cols)


#This will make our bray cutris distance matrix for beta diversity calculations.
bray.dist.matrix <- as.data.frame(as.matrix(avgdist(q1.table, sample = 1000, meanfun = median, transf = sqrt, iterations = 999))) #default is Bray-Curtis calculated by vegdist()
write.table(x = bray.dist.matrix, file = "median_bray_curtis_sqrt_q1_single_rarefied_to_1000.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


















