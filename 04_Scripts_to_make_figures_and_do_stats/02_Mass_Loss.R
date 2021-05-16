# Combine datasets for mass loss data and then analyze it!
# Kendra E Walters
# August 14th 2020

### 
require(data.table)
require(ggplot2)
require(car)


##### ____________ COMBINE THE DATASETS ___________________ ######################


setwd("/Users/walters_kendra/Google Drive/Dispersal_Ch_3/02_Data/05_Mass_Loss/")


# Convert envelope weights to their dry weight (they lose about 5% of mass during drying process)
envelope.weights <- as.data.frame(fread("All_empty_envelope_weights - Sheet1.csv"))
envelope.dry.wet.ratio <- as.data.frame(fread("T2_envelope_dry_wet_ratio - Sheet1.csv"))
envelope.weights$dry_empty_envelope_g <- envelope.weights$`Empty Envelope (g)` * mean(envelope.dry.wet.ratio$`Dry:Wet Ratio`)

# Read in meta data and initial litterbag weights
initial.litter.weights <- as.data.frame(fread("Initial_Litter_Weights - Sheet1.csv"))
metadata <- as.data.frame(fread("../01_Master_Metadata - Sheet1.tsv"))
initial.meta <- merge(initial.litter.weights, metadata, by = "Old_Sample_ID", all.x = TRUE)

# Pare down initial litterbag dataset
keep <- c("Old_Sample_ID", "New_Sample_ID", "Initial Green Litter Mass (g)", "Substrate", "Treatment", "Time_Point", "Transect", "Ecosystem", "Location")
initial.meta <- initial.meta[ , keep]
initial.meta <- initial.meta[!(is.na(initial.meta$`Initial Green Litter Mass (g)`)), ]

# Getting T0 dry:wet ratio for grass and shrub
t0.dry.wet <- as.data.frame(fread("T0_Dry_Wet_Ratio - Sheet1.csv")) # this mass is of JUST the litter so you don't need to convert envelope weights or anything
grass.t0.dry.wet <- median(t0.dry.wet[t0.dry.wet$`Litter Type` == "Grass", ]$Dry_to_Wet_Ratio)
shrub.t0.dry.wet <- median(t0.dry.wet[t0.dry.wet$`Litter Type` == "Shrub", ]$Dry_to_Wet_Ratio)

# Read in all timepoints for Person 1 and combine
t1.person.1 <- as.data.frame(fread("T1_Person_1 - Sheet1.csv"))
t2.person.1 <- as.data.frame(fread("T2_Person_1 - Sheet1.csv"))
t4.person.1 <- as.data.frame(fread("T4_Person_1 - Sheet1.csv"))

person.1 <- rbind(t1.person.1, t2.person.1)
person.1 <- rbind(person.1, t4.person.1)
names(person.1)[5] <- "Initials_Person_1"

# Read in all timepoints for Person 2 and combine
t1.person.2 <- as.data.frame(fread("T1_Person_2 - Sheet1.csv"))
t2.person.2 <- as.data.frame(fread("T2_Person_2 - Sheet1.csv"))
t4.person.2 <- as.data.frame(fread("T4_Person_2 - Sheet1.csv"))

person.2 <- rbind(t1.person.2, t2.person.2)
person.2 <- rbind(person.2, t4.person.2)
names(person.2)[4] <- "Initials_Person_2"

# Read in all timepoints for dry envelope + litter weights and combine
t1.dry.env.and.litter <- as.data.frame(fread("T1_dry_envelope_and_litter_weight - Sheet1.csv"))
t2.dry.env.and.litter <- as.data.frame(fread("T2_dry_envelope_and_litter_weight - Sheet1.csv"))
t4.dry.env.and.litter <- as.data.frame(fread("T4_dry_envelope_and_litter_weight - Sheet1.csv"))

dry.env.and.litter <- rbind(t1.dry.env.and.litter, t2.dry.env.and.litter)
dry.env.and.litter <- rbind(dry.env.and.litter, t4.dry.env.and.litter)

# calculate dry:wet ratio
combo <- merge(dry.env.and.litter, envelope.weights, by = "SampleID", all.x = TRUE)
combo.combo <- merge(combo, person.1, by = "SampleID")
combo.combo.combo <- merge(combo.combo, person.2, by = "SampleID")

combo.combo.combo$Wet_Litter_In_Envelope <- combo.combo.combo$`Envelope + Wet Litter (g)` - combo.combo.combo$`Empty Envelope (g)`
combo.combo.combo$Dry_Litter_In_Envelope <- combo.combo.combo$`Dry Envelope + Litter (g)` - combo.combo.combo$dry_empty_envelope_g

combo.combo.combo$Dry_Wet_Ratio <- combo.combo.combo$Dry_Litter_In_Envelope / combo.combo.combo$Wet_Litter_In_Envelope

# calculate final dry litter weight
combo.combo.combo$Wet_Litter_All <- combo.combo.combo$`Full Litterbag (g)` - combo.combo.combo$`Empty Litterbag (g)`
combo.combo.combo$Dry_Litter_All <- combo.combo.combo$Wet_Litter_All * combo.combo.combo$Dry_Wet_Ratio

# calculate initial dry litter weight
initial.meta$Initial_Dry_Weight <- ifelse(initial.meta$Substrate == "Grass Litter", 
                                          initial.meta$`Initial Green Litter Mass (g)` * grass.t0.dry.wet, 
                                          initial.meta$`Initial Green Litter Mass (g)` * shrub.t0.dry.wet) 


# merge initial with final measurements
final <- merge(combo.combo.combo, initial.meta, by.x = "SampleID", by.y = "New_Sample_ID")

# calculate mass loss
final$Mass_Loss <- 1 - final$Dry_Litter_All / final$Initial_Dry_Weight

final$Time_Point <- as.factor(final$Time_Point)

pdf("../../04_Figures/01_Figure_Exploration/02_mass_loss_all.pdf")
ggplot(data = final) + geom_boxplot(aes(x = paste(Substrate, Time_Point), y = Mass_Loss, fill = Treatment)) + 
  theme_classic() + xlab("Litter Substrate & Time Point") + ylab("Mass Loss")
dev.off()
final$Treatment <- as.factor(final$Treatment)
final$Substrate <- as.factor(final$Substrate)
final$Time_Point <- as.factor(final$Time_Point)
mod <- lm(formula = Mass_Loss ~ Substrate + Treatment, data = final)
Anova(mod, type = 3)

ggplot(final) + geom_point(aes(x = Initial_Dry_Weight, y = Dry_Litter_All, color = Substrate))
+ theme_classic() + geom_abline(slope = 1)



