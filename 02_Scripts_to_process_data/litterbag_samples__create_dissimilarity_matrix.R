# Script to make a BC dissimilarity matrix given output from QIIME
# Dec. 23rd 2020 
# Kendra E Walters

# Room of ()ments
require(vegan)
require(data.table)

# Get ready for working with the data
setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_GitHub_Data_and_Scripts/01_Raw_data/06_Community_Composition/06_Output_from_qiime/Q2_single_16S/")

input.data <- "12_table_filtered.tsv"
output.data <- "/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_GitHub_Data_and_Scripts/03_Processed_data/litterbag_samples__16S_single__BC_median_r1000.tsv"
output.rda <- "/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_GitHub_Data_and_Scripts/03_Processed_data/litterbag_samples__16S_single__BC_median_r1000.rda"

# Read in the data
q1.table <- as.data.frame(fread(input.data)) 
row.names(q1.table) <- q1.table$`#OTU ID` # make otuID the row names
q1.table <- q1.table[ , !(names(q1.table) %in% c("#OTU ID"))]
q1.table <- as.data.frame(t(q1.table)) # vegan expects rows = samples and columns = species

# Choose a rarefaction depth
sort(rowSums(q1.table)) # choose a rarefaction depth
hist(rowSums(q1.table), breaks = 80) # oof, the sequencing did not go well here
nrow(q1.table[rowSums(q1.table) < 1000, ]) / nrow(q1.table) # see how many samples you might lose

# Calculate matrix
X <- 1000 # rarefaction depth
bray.dist <- avgdist(q1.table, sample = X, meanfun = median, transf = sqrt, iterations = 999) #default is Bray-Curtis calculated by vegdist()
bray.dist.df <- as.data.frame(as.matrix(bray.dist))
write.table(x = bray.dist.df, output.data, sep = "\t", row.names = TRUE)

save(bray.dist, file = output.rda)


















