# Whyyyy does my q2 sequencing suck so much??
# January 29th 2021
# Kendra E Walters


# Room of ()ments
require(data.table)
require(stringr)
require(tidyverse)

##### ________ This part is to look at the OTU table pre- and post-filtering in QIIME to see what it is taking out ________ ########

# Reading in and creating the right data format
setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_Data/06_Community_Composition/06_Output_from_qiime/Q2_single_16S/")

QIIME.pre <- fread("14_table_not_filtered.tsv", data.table = FALSE)
taxon.pre <- fread("15_taxonomy_not_filtered.tsv", data.table = FALSE)

QIIME.pre <- merge(QIIME.pre, taxon.pre, by.x = "#OTU ID", by.y = "Feature ID", all.x = TRUE, all.y = FALSE)

QIIME.post <- fread("14_table_filtered_with_taxonomy_filtered.tsv", data.table = FALSE)

# Looking at the difference between the two dataframes
taxon.pre <- QIIME.pre$Taxon
taxon.post <- QIIME.post$taxonomy

QIIIME.taxon.filtered <- setdiff(taxon.pre, taxon.post) # the UNIQUE taxon IDs from the 1772 OTU IDs that were filtered out

# So looks like it is filtering out: 
  # "k__Bacteria"
  # "Unassigned"
  # "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__mitochondria; g__; s__"
  # "k__Bacteria; p__Cyanobacteria; c__Chloroplast (with many lower classifications)


##### ________ Here we are, filtering in R, and exploring chlorlplast and phylum assignment frequencies _____######
# Loading and merging data
setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_Data/06_Community_Composition/06_Output_from_qiime/Q2_single_16S/")

table <- fread("14_table_not_filtered.tsv", data.table = FALSE)
taxa <- fread("15_taxonomy_not_filtered.tsv", data.table = FALSE)

table.taxa <- merge(table, taxa, by.x = "#OTU ID", by.y = "Feature ID", all.x = TRUE, all.y = FALSE)
table.taxa$`#OTU ID` <- NULL

# Row sum to give # of sequences 
sumsum <- rowSums(table.taxa[ , !(names(table.taxa) %in% c("Taxon", "Confidence"))])
sum.reads.taxa <- data.frame("Total_Read_Count" = sumsum, "Taxon" = table.taxa$Taxon)


# What percent of reads are identified to phylum? 
sum(sum.reads.taxa[ grepl("p__[A-Za-z]", sum.reads.taxa$Taxon), ]$Total_Read_Count) / sum(sum.reads.taxa$Total_Read_Count) # 83.2%  of reads identified to phylum

# What percent are classified as chloroplasts? 
sum(sum.reads.taxa[ grepl("Chloroplast", sum.reads.taxa$Taxon), ]$Total_Read_Count) / sum(sum.reads.taxa$Total_Read_Count) # 63% of reads are chloroplasts
sum(sum.reads.taxa[ grepl("mitochondria", sum.reads.taxa$Taxon), ]$Total_Read_Count) / sum(sum.reads.taxa$Total_Read_Count) # 15.1% are mitochondria

# of the non-Chloroplast data, how many are identified to phyllum
no.chloro.mito <- sum.reads.taxa %>% filter(!(str_detect(Taxon, "Chloroplast|mitochondria")))
sum(no.chloro.mito[grepl("p__[/[A-Za-z]", no.chloro.mito$Taxon), ]$Total_Read_Count) / sum(no.chloro.mito$Total_Read_Count) # 23.4%  of reads identified to phylum


# Looking at reads/sample without chloroplast DNA
table.taxa.cleaned <- table.taxa %>% filter(!(str_detect(Taxon, "Chloroplast|mitochondria")) & str_detect(Taxon, "p__[\\[A-Za-z]"))
#[!(grepl("Chloroplast", table.taxa$Taxon)), ] %>% filter(grepl("p__[A-Za-z]", Taxon))

cleaned.sample.read.counts <- sort(colSums(table.taxa.cleaned[ , !(names(table.taxa.cleaned) %in% c("Taxon", "Confidence"))])) %>% as.data.frame
names(cleaned.sample.read.counts) <- c("Sequence_Count_by_R")
length(cleaned.sample.read.counts[cleaned.sample.read.counts > 600]) # 48 ---- seems comparable to QIIME
48/317 # 15% .... maybe worse than QIIME??


# Read in QIIME cleaned sequences / sample data
qiime.cleaned.sample.read.counts <- fread("sample-frequency-detail-from-QIIME2-View.csv", data.table = FALSE)
names(qiime.cleaned.sample.read.counts) <- c("SampleID", "Sequence_Count_by_QIIME")

cleaned.sample.read.counts.together <- merge(cleaned.sample.read.counts, qiime.cleaned.sample.read.counts, 
                                             by.x = "row.names", by.y = "SampleID")
cleaned.sample.read.counts.together$Difference <- cleaned.sample.read.counts.together$Sequence_Count_by_R - cleaned.sample.read.counts.together$Sequence_Count_by_QIIME
sum(cleaned.sample.read.counts.together$Difference) ## WooHOO!! They are exactly the same





