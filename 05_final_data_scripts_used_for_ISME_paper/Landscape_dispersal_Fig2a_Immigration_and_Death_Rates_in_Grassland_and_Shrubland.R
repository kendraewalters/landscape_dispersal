## Figure 2A Immigration + Death Rates in Grassland and Shrubland
## Created on July 19 2021 by Kendra E Walters
## Edited by Kristin Barbour with updated flow cytometry data (only subtracted negative controls to closed bags)
## Also updated the immigration rate (accumulation) plot to just include a linear model rather than the dynamic model


# Room of ()ments
require(ggplot2)
require(vegan)
library(plyr)
require(tidyverse)
require(conflicted)
require(egg)
library(patchwork)
library(car)

conflict_prefer(name = "select", winner = "dplyr")
conflict_prefer(name = "filter", winner = "dplyr")
conflict_prefer(name = "mutate", winner = "dplyr")
conflict_prefer(winner = "dplyr", name = "summarize")
conflict_prefer(winner = "dplyr", name = "summarise")

## set working directory
setwd("~/kbarbou1@uci.edu - Google Drive/My Drive/Professional_development/Kendra_Ch3_Manuscript/03_Analysis/")

###### SETTING UP OUR DATA ######

# Load the data
# load("Kendra_processed_data/flow_cytometry_cleaned.rda")
load("Kendra_processed_data/metadata_cleaned.rda")

## update as of 01/31/2025, altered the way negative controls were subtracted from all samples. See "clean_flow_cytometry_data_KB.R" for more information
cleaned.cell.counts.KB <- read.csv("cleaned_flow_data_KB.csv", sep=",")

# Set up a key to change time point into # of days
days.key <- c("0" = 0, "1" = 12, "2" = 26, "3" = 47, "4" = 54)

# Clean up our data to just be what we need for this analysis 
cell.counts.use <- cleaned.cell.counts.KB %>% 
  left_join(metadata.together %>% 
              mutate(`Random_._Flow` = as.character(`Random_#_Flow_Cytometry`))) %>% # combine with metadata using Random_#_Flow
  mutate(Treatment = ifelse(is.na(Treatment), "Death", Treatment)) %>% # add metadata for T0 samples that weren't included in the overall metadata file
  mutate(Time_Point = ifelse(is.na(Time_Point), 0, Time_Point)) %>% 
  mutate(Ecosystem = ifelse(grepl("g", `Random_._Flow`), "Grassland",
                            ifelse(grepl("s", `Random_._Flow`), "Shrubland", 
                                   ifelse(grepl("n", `Random_._Flow`), "Negative", Ecosystem)))) %>%
  select(number_cells_on_glass_slide_cleaned, Treatment, Ecosystem, Time_Point) %>% # reduce to only columns we need
  filter(Ecosystem != "Negative") %>% # remove this negative control (should have been removed prior to this script)
  mutate(Time_as_days = days.key[as.character(Time_Point)]) %>% # make column that lists # of days
  mutate(cells_per_sq_cm = number_cells_on_glass_slide_cleaned / (7.5 * 2.5))

# Separate into our two ecosystems
grass.cells <- cell.counts.use %>% filter(Ecosystem == "Grassland")
shrub.cells <- cell.counts.use %>% filter(Ecosystem == "Shrubland")
all.cells <- cell.counts.use

## just do a quick summary of the abundances
cell.summary <- all.cells %>%
  as.data.frame() %>% 
  # filter(Treatment=="Death") %>% 
  ddply(., c("Treatment", "Ecosystem", "Time_Point"), summarise, cells = mean(cells_per_sq_cm))

## Does overall abundance on the glass slides change overtime?
abundance.time.model <- all.cells %>%
  dplyr::filter(Treatment == "Open") %>%
  lm(cells_per_sq_cm ~ Ecosystem*Time_as_days, data = .)

Anova(abundance.time.model, type = "II") ## Yes, timepoint is significant

###### DEATH RATES ######

# Calculating death rate for the grassland death slides
intercept.grass <- grass.cells %>% 
  filter(Treatment == "Death", 
         Time_as_days == 0) %>% 
  summarize(mean = mean(log(cells_per_sq_cm))) %>% pull(mean)

r.grass <- lm(I(log(cells_per_sq_cm) - intercept.grass) ~ Time_as_days + 0, 
              grass.cells %>% filter(Treatment == "Death")) %>% coef*-1 # This is calculating death rate for the grassland death slides

## extract standard error
r.grass.std.error <- summary(lm(I(log(cells_per_sq_cm) - intercept.grass) ~ Time_as_days + 0, 
                        grass.cells %>% filter(Treatment == "Death"))) %>% 
  .[["coefficients"]] %>% 
  .[,2] #0.002938828

# Calculating death rate for the shrub death slides
intercept.shrub <- shrub.cells %>% 
  filter(Treatment == "Death", 
         Time_as_days == 0) %>% 
  summarize(mean = mean(log(cells_per_sq_cm))) %>% pull(mean)

r.shrub <- lm(I(log(cells_per_sq_cm) - intercept.shrub) ~ Time_as_days + 0, 
              shrub.cells %>% filter(Treatment == "Death")) %>% coef*-1 # this is calculating death rate for the shrub death slides

## extract standard error
r.shrub.std.error <- summary(lm(I(log(cells_per_sq_cm) - intercept.shrub) ~ Time_as_days + 0, 
                                shrub.cells %>% filter(Treatment == "Death"))) %>% 
  .[["coefficients"]] %>% 
  .[,2] #0.00236632

# creating a mini data frame to give our parameters for the linear models (since we are picky about our intercepts here)
lm.parameters.death <- data.frame("Intercept" = c(intercept.grass, intercept.shrub), 
                                  "Slope" = c(r.grass[[1]], r.shrub[[1]]), 
                                  "Ecosystem" = c("Grassland", "Shrubland"))

# plotting our death rate data
A <- ggplot(data = cell.counts.use %>% filter(Treatment == "Death") %>% 
              filter(Ecosystem!="Negative") %>% 
         mutate(Time_as_days= as.integer(Time_as_days)), 
       aes(x = as.integer(Time_as_days), y = log(cells_per_sq_cm), 
           color = Ecosystem)) +
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_abline(data = lm.parameters.death, aes(intercept = Intercept, slope = -Slope, color = Ecosystem), linewidth = 1) + 
  scale_color_manual(values = c("Grassland" = "#508104", "Shrubland" = "#db8200") 
                    ) + 
  facet_wrap(~ Treatment, scales = "free_y") +
  theme_test() + 
  labs(y = bquote("Bacterial cells / " ~cm^2), x = "Number of Days", color = "Ecosystem") + 
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=17, face = "bold", 
                                  margin = margin(4, 0, 4, 0)),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour="black", size = 12), 
        axis.text.y = element_text(colour="black", size = 12), 
        axis.title.x = element_text(vjust = 0.35, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14), 
        legend.position = "none") +
  scale_x_continuous(limits = c(0,55), n.breaks = 6)

plot(A)

### Does death rate vary between ecosystems? No.
death.ancova <- cell.counts.use %>% filter(Treatment == "Death") %>% 
  filter(Ecosystem!="Negative") %>% 
  lm(formula = log(cells_per_sq_cm) ~ Ecosystem*Time_as_days, data = .) %>% 
  Anova(type = "II")

death.ancova ## Time x Ecosystem P = 0.167

###### IMMIGRATION RATES ######

# Calculating the death rate and immigration rate coefficients WITH measured death rates
(nls.shrub <- nls(cells_per_sq_cm ~  
                    Im/r.shrub * (1 - exp(-r.shrub * Time_as_days)), 
                  data = shrub.cells[shrub.cells$Treatment == "Open", ],
                  start = list(Im = 2e5)) %>% summary) ## Im = 953.07, std error = 79.15

(nls.grass <- nls(cells_per_sq_cm ~ 
                    Im/r.grass * (1 - exp(-r.grass * Time_as_days)), 
                  data = grass.cells[grass.cells$Treatment == "Open", ],
                  start = list(Im = 2e5)) %>% summary) ## Im = 1178.3, std error = 168.8

## does immigration rate vary between ecosystems? Using an ANCOVA (assumes a linear model but should be fine)
intercept.grass.im <- grass.cells %>% 
  filter(Treatment == "Open",
         Time_Point == 1) %>% 
  summarize(mean = mean(log(cells_per_sq_cm))) %>% pull(mean)


intercept.shrub.im <- shrub.cells %>% 
  filter(Treatment == "Open", 
         Time_Point == 1) %>% 
  summarize(mean = mean(log(cells_per_sq_cm))) %>% pull(mean)

immigration.ancova <- cell.counts.use %>% dplyr::filter(Treatment == "Open") %>% 
  mutate(log_abundance_minus_intercept = ifelse(Ecosystem == "Shrubland",
                                                log(cells_per_sq_cm) - intercept.shrub.im,
                                                ifelse(Ecosystem == "Grassland",
                                                       log(cells_per_sq_cm) - intercept.grass.im, NA))) %>% 
  lm(formula = log_abundance_minus_intercept ~ Ecosystem*Time_as_days, data = .) %>% 
  Anova(type = "II") # Not significant!!

immigration.ancova ## Time x Ecosystem P = 0.4691

## Calculate average immigration rate using with the average death rate across ecosystems (because neither immigration rate nor death rate differs between ecosystem)
## calculate average death rate
r.all <- (r.grass + r.shrub)/2

## immigration rate across all samples using the average death rate
(nls.all <- nls(cells_per_sq_cm ~
                  Im/r.all * (1 - exp(-r.all * Time_as_days)),
                data = all.cells[all.cells$Treatment == "Open", ],
                start = list(Im = 2e5)) %>% summary) ## Im = 1,060.57, std error = 89.73 

## calculate linear model to have accumulation rates
## need to first add a T0 = 0 for the open samples, need to include that all glass slides started out as sterile with no cells on them
grass.cells.accumulation <- add_row(grass.cells, Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>%
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Grassland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0)
  
shrub.cells.accumulation <- add_row(shrub.cells, Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0) %>% 
  add_row(., Ecosystem = "Shrubland", Treatment = "Open", number_cells_on_glass_slide_cleaned = 0, Time_Point = 0, Time_as_days = 0, cells_per_sq_cm = 0)


## linear model, forcing the model through 0 intercept (will just be calculating the slope)
accumulation.reg.grass <- grass.cells.accumulation %>% 
  filter(Treatment=="Open") %>% 
  lm(formula = cells_per_sq_cm ~ 0 + Time_as_days, data = .)

accumulation.grass.coeff <- coefficients(accumulation.reg.grass)
accumulation.grass.intercept <- 0
accumulation.grass.slope <- accumulation.grass.coeff[1]

## linear model, forcing the model through 0 intercept (will just be calculating the slope)
accumulation.reg.shrub <- shrub.cells.accumulation %>% 
  filter(Treatment=="Open") %>% 
  lm(formula = cells_per_sq_cm ~ 0 + Time_as_days, data = .)

accumulation.shrub.coeff <- coefficients(accumulation.reg.shrub)
accumulation.shrub.intercept <- 0
accumulation.shrub.slope <- accumulation.shrub.coeff[1]

## test whether accumulation rates vary between ecosystems
accumulation.ancova <- rbind(grass.cells.accumulation, shrub.cells.accumulation) %>% 
  filter(Treatment == "Open") %>% 
  lm(formula = cells_per_sq_cm ~ Time_as_days*Ecosystem, data = .) %>% 
  Anova(type = "II")

accumulation.ancova

## combine datasets for figure
both.cells.disp <- grass.cells.accumulation %>% 
  mutate(Ecosystem = c("Grassland")) %>% 
  bind_rows(shrub.cells.accumulation %>% 
              mutate(Ecosystem = c("Shrubland"))) %>% 
  filter(Treatment == "Open") %>% 
  mutate(Days = as.integer(Time_as_days))

# Graph the new curves that we made
cols <- c("Grassland" = "#508104", "Shrubland"  = "#db8200") # slightly darker set of orange + green

## decided to just include the accumulation curves
B <- ggplot(data = both.cells.disp, 
            aes(x = Days, y = cells_per_sq_cm, color = Ecosystem)) +
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_abline(intercept = accumulation.grass.intercept, 
              slope = accumulation.grass.slope, 
              color = cols[["Grassland"]], 
              linewidth = 1,
              linetype = "solid") +
  geom_abline(intercept = accumulation.shrub.intercept, 
              slope = accumulation.shrub.slope, 
              color = cols[["Shrubland"]], 
              linewidth = 1,
              linetype = "solid") +
  scale_color_manual(values = cols) + 
  facet_wrap(~ Treatment, scales = "free_y") +
  scale_x_continuous(limits = c(0,55), n.breaks = 6) +
  theme_test() + 
  labs(y = bquote("log(Bacterial cells / " ~cm^2), x = "Number of Days", color = "Ecosystem") + 
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=17, face = "bold", 
                                  margin = margin(4, 0, 4, 0)),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour="black", size = 12), 
        axis.text.y = element_text(colour="black", size = 12), 
        axis.title.x = element_text(vjust = 0.35, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14), 
        legend.title = element_text(size= 16), 
        legend.text = element_text(size = 15))

plot(B)

### EXPORT GRAPH
B + A + plot_layout(guides = "collect")
ggsave(plot=last_plot(), file="Immigration_Death_Rates_in_Grassland_and_Shrubland_v5.pdf", width = 8.79, height = 3.39, units = "in")

###### Comparison to litter ######
all_immigration_rate <- nls.all$parameters[[1]]

## landscape scale (don't split by ecosystem)
read.delim("Kendra_processed_data/BacterialAbundanceLomaRidgeLitterbagsFromClaudia.txt") %>% 
  filter(timepoint!=0,
         moTreatment == "ambient") %>% 
  mutate(`BACavgcount.g.dry.weight` = as.numeric(`BACavgcount.g.dry.weight`), 
         AbundancePerBag = `DryMass..g.` * `BACavgcount.g.dry.weight`, 
         AbundancePer_cm2 = AbundancePerBag / (13 * 13)) %>% # assuming a 15x15cm bag with 1 cm sealed with duct-tape on each of the 4sides (results don't change between 13x13 or 13x14 - all 0.04%) 
  # group_by(LitterOrigin) %>% 
  summarize(mean = mean(AbundancePer_cm2)) %>% 
  mutate(immigration_rate = all_immigration_rate,
    immigration_percent = (all_immigration_rate / mean) *100) %>% 
  print(.)
