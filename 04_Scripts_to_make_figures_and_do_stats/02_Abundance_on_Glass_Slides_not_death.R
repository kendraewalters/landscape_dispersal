# Abundance measurements on slides, not death (just dispersal)
# 7/2/20
# Kendra E Walters

# Room of ()ments
require(data.table)
require(ggplot2)


# Read in the data
setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_Data/04_Flow_Cytometry/")
flow <- as.data.frame(fread("adjusted_gates__dispersal2__T0_T1_T2_T3_T4__just_needed_samples.csv"))
meta <- as.data.frame(fread("./../01_Master_Metadata - Sheet1.tsv"))


# Match the SampleID to the Random_# in the flow dataset
flow.meta <- merge(flow, meta, by.x = "Random_#", by.y = "Random_#_Flow", all.x = TRUE, all.y = FALSE)

# Add timepoints for the negatives (in lab controls, not in field controls) and T0 death slides:
flow.meta$Time_Point <- ifelse(grepl("T1", flow.meta$Sample), "1", 
                               ifelse(grepl("T2", flow.meta$Sample), "2", 
                                      ifelse(grepl("T3", flow.meta$Sample), "3", 
                                             ifelse(grepl("T4", flow.meta$Sample), "4", 
                                                    ifelse(flow.meta$Specimen == "T0__Death__61920", "0", flow.meta$Time_Point)))))

flow.meta$Treatment <- ifelse(grepl("T0|T1|T2|T3|T4", flow.meta$Sample), "Lab Negative",
                              ifelse(flow.meta$Specimen == "T0__Death__61920", "Death",flow.meta$Treatment))

flow.meta$Substrate <- ifelse(flow.meta$Specimen == "T0__Death__61920" & grepl("g", flow.meta$Sample), "Cells from Grass", 
                              ifelse(flow.meta$Specimen == "T0__Death__61920" & grepl("s", flow.meta$Sample), "Cells from Shrub", flow.meta$Substrate))

flow.meta$Ecosystem <- ifelse(flow.meta$Specimen == "T0__Death__61920" & grepl("g", flow.meta$Sample), "Grassland", 
                             ifelse(flow.meta$Specimen == "T0__Death__61920" & grepl("s", flow.meta$Sample), "Shrubland", flow.meta$Ecosystem))


# Calculating the number of cells based on dilution and volumes of liquids added
flow.meta$Cell_Calc <- flow.meta$`E4 Abs. Count` * flow.meta$Dilution * flow.meta$End_Volume_uL * (1/flow.meta$Mid_Volume_uL) * flow.meta$Starting_Volume_uL

# Getting an average for negatives (lab and field) for each time point and visualize it
time.points <- c("1", "2", "3", "4")

for (i in 1:length(time.points)) {
  subset <- flow.meta[flow.meta$Time_Point == time.points[[i]], ]
  lab_control <- mean(subset[subset$Treatment == "Lab Negative", ]$Cell_Calc)
  field_control <- mean(subset[subset$Treatment == "Closed", ]$Cell_Calc)
  field_contamination <- field_control - lab_control
  print(paste0("For timepoint ", time.points[[i]], ", your lab contamination was ", lab_control, ", your field contamination was ", field_control, ", and the difference is ", field_contamination, "."))
}

pdf("./../../04_Figures/01_Figure_Exploration/01_first_flow_boxplot_to_look_at_contamination_levels.pdf")
ggplot(data = flow.meta) + geom_boxplot(aes(x = Treatment, y = Cell_Calc, fill = Time_Point)) + theme_classic()
dev.off()

# Contamination looks GREAT (by that, I mean low) 
# Let's subtract the field control averages from each samples for each time point

t0.neg <- mean(flow.meta[flow.meta$Time_Point == "0" & flow.meta$Treatment == "Closed", ]$Cell_Calc) # this one is lab negative cause we don't have field negative for T0
t1.neg <- mean(flow.meta[flow.meta$Time_Point == "1" & flow.meta$Treatment == "Closed", ]$Cell_Calc)
t2.neg <- mean(flow.meta[flow.meta$Time_Point == "2" & flow.meta$Treatment == "Closed", ]$Cell_Calc)
t3.neg <- mean(flow.meta[flow.meta$Time_Point == "3" & flow.meta$Treatment == "Closed", ]$Cell_Calc)
t4.neg <- mean(flow.meta[flow.meta$Time_Point == "4" & flow.meta$Treatment == "Closed", ]$Cell_Calc)

for.analysis <- flow.meta[flow.meta$Treatment != "Closed" & flow.meta$Treatment != "Lab Negative", ]

for.analysis$Adjusted_Cell_Calc <- ifelse(for.analysis$Time_Point == "1", for.analysis$Cell_Calc - t1.neg, 
                                       ifelse(for.analysis$Time_Point == "2", for.analysis$Cell_Calc - t2.neg, 
                                              ifelse(for.analysis$Time_Point == "3", for.analysis$Cell_Calc - t3.neg,
                                                     ifelse(for.analysis$Time_Point == "4", for.analysis$Cell_Calc - t4.neg, 
                                                            ifelse(for.analysis$Time_Point == "0", for.analysis$Cell_Calc - t1.neg, "YIKES! Try again...")))))


# For the analysis, you will want the timepoint as weeks
for.analysis$Weeks <- ifelse(for.analysis$Time_Point == 0, 0, 
                             ifelse(for.analysis$Time_Point == 1, 2,
                                    ifelse(for.analysis$Time_Point == 2, 4,
                                           ifelse(for.analysis$Time_Point == 3, 7,
                                                  ifelse(for.analysis$Time_Point == 4, 8, "eeeek... check your script again.")))))
for.analysis$Weeks <- as.numeric(for.analysis$Weeks)


# Let's first analyze the death slides
death <- for.analysis[for.analysis$Treatment == "Death", ]
death$Time_Point <- as.numeric(death$Time_Point)
death$Adjusted_Cell_Calc <- as.numeric(death$Adjusted_Cell_Calc)

summary(lm(formula = Adjusted_Cell_Calc ~ Time_Point, data = death))
summary(lm(formula = log(Adjusted_Cell_Calc) ~ Time_Point, data = death))


summary(lm(formula = log10(Adjusted_Cell_Calc) ~ Weeks, data = death)) # -0.07673


intercept <- mean(log10(death[death$Weeks == 0, ]$Adjusted_Cell_Calc))

lm_death <- lm(I(log10(Adjusted_Cell_Calc) - intercept) ~ Weeks + 0, death) # estimate: -0.08318
summary(lm_death)

1 - 10^-0.08318 # 0.1743043 === 17.43% per week


ggplot() + 
  geom_point(data = death, aes(x = `Weeks`, y = log10(Adjusted_Cell_Calc), color = Ecosystem)) + 
  theme_classic() +
  xlab("Weeks") + 
  ylab("log10(Abundance)") + 
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
  geom_abline(intercept = intercept, slope = -0.08318, linetype = 2)


# Now let's analyze the dispersal slides
for.analysis$Adjusted_Cell_Calc <- as.numeric(for.analysis$Adjusted_Cell_Calc)
for.analysis.open <- for.analysis[for.analysis$Treatment == "Open", ]

pdf("./../../04_Figures/01_Figure_Exploration/03_abundance_by_ecosystem_boxplot.pdf")
ggplot(data = for.analysis.open) + geom_boxplot(aes(x = Ecosystem, y = Adjusted_Cell_Calc, fill = Ecosystem)) + 
  theme_classic()
dev.off()

pdf("./../../04_Figures/01_Figure_Exploration/04_abundance_by_ecosystem_and_timepoint_boxplot.pdf")
ggplot(data = for.analysis.open) + geom_boxplot(aes(x = Ecosystem, y = Adjusted_Cell_Calc, fill = Time_Point)) + 
  theme_classic()
dev.off()

ggplot(data = for.analysis.open) + geom_boxplot(aes(x = Time_Point, y = Adjusted_Cell_Calc, fill = Time_Point)) + 
  theme_classic()

pdf("./../../04_Figures/01_Figure_Exploration/05_abundance_by_transect_boxplot.pdf")
ggplot(data = for.analysis.open) + geom_boxplot(aes(x = Transect, y = Adjusted_Cell_Calc, fill = Transect)) + 
  theme_classic()
dev.off()

pdf("./../../04_Figures/01_Figure_Exploration/06_abundance_by_transect_and_timepoint_boxplot.pdf")
ggplot(data = for.analysis.open) + geom_boxplot(aes(x = Transect, y = Adjusted_Cell_Calc, fill = Time_Point)) + 
  theme_classic()
dev.off()

pdf("./../../04_Figures/01_Figure_Exploration/06_grassland_abundance_by_transect_and_timepoint_boxplot.pdf")
ggplot(data = for.analysis.open[for.analysis.open$Ecosystem == "Grassland", ]) + geom_boxplot(aes(x = Transect, y = Adjusted_Cell_Calc, fill = Time_Point)) + 
  theme_classic()
dev.off()

pdf("./../../04_Figures/01_Figure_Exploration/06_shrubland_abundance_by_transect_and_timepoint_boxplot.pdf")
ggplot(data = for.analysis.open[for.analysis.open$Ecosystem == "Shrubland", ]) + geom_boxplot(aes(x = Transect, y = Adjusted_Cell_Calc, fill = Time_Point)) + 
  theme_classic()
dev.off()

pdf("./../../04_Figures/01_Figure_Exploration/06_abundance_by_transect_and_ecosystem_boxplot.pdf")
ggplot(data = for.analysis.open) + geom_boxplot(aes(x = Transect, y = Adjusted_Cell_Calc, fill = Ecosystem)) + 
  theme_classic()
dev.off()
