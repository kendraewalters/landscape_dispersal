# Analyzing death slides for Chapter 3, Dispersal Experiment
# Kendra E Walters
# June 23rd, 2020

# Room of ()ments
require(data.table)
require(ggplot2)


# Read in the data
setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_Data/Flow_Cytometry/")
death <- as.data.frame(fread("all_death_slides__cleaned.csv"))


more.death <- flow.meta.cleaned[flow.meta.cleaned$Treatment == "Death", ]


# Create our useful columns
death$litter_source <- ifelse(grepl(pattern = "g", death$Sample), "Grassland", 
                              ifelse(grepl("s", death$Sample), "Shrubland", "whoooooops!"))
death$Avg_Reading_on_Flow <- rowMeans(death[ , c("E4 Abs. Count", "E5 Abs. Count", "E8 Abs. Count")])

death$Cells_on_Slide <- death$Avg_Reading_on_Flow * death$Dilution_Factor * death$End_Volume_uL * (1/death$Mid_Volume_uL) * death$Starting_Volume_uL

# Looking at averages
aggregate(x = death$Cells_on_Slide, list(death$litter_source), mean)

# TEMP working space:
death.keep <- death[, c("litter_source", "Cells_on_Slide")]
death.keep$Time_Point <- c(rep("0", nrow(death.keep)))
names(death.keep) <- c("Ecosystem", "Cells_on_Slide", "Time_Point")

a <- death.keep
a$Transect <- c(rep("A", nrow(death.keep)))
b <- death.keep
b$Transect <- c(rep("B", nrow(death.keep)))
c <- death.keep
c$Transect <- c(rep("C", nrow(death.keep)))

to <- rbind(a, b)
tog <- rbind(to, c)

more.death.keep <- more.death [ , c( "Ecosystem", "Cells_on_Slide_minus_LO", "Time_Point", "Transect")]
names(more.death.keep) <- c("Ecosystem", "Cells_on_Slide", "Time_Point", "Transect")
more.death.keep.minus <- more.death.keep
more.death.keep.minus$Transect <- NULL
simple.death.all <- rbind(death.keep, more.death.keep.minus)

ggplot(data = simple.death.all) + geom_boxplot(aes(x = Ecosystem, y = Cells_on_Slide, fill = Time_Point)) + 
  theme_classic()

aggregate(simple.death.all$Cells_on_Slide, by = list(simple.death.all$Ecosystem, simple.death.all$Time_Point), mean)
631237.33/1105130.48 # 57.12 % died in grassland (= 81.5 % per month, compared to 60% last time)
78098.95/144782.10 # 53.94 % died in shrubland (= 78.84 % per month)

# Plotting the data
ggplot(data = death) + geom_boxplot(aes(x = litter_source, y = Cells_on_Slide)) + 
  theme_classic()
