# make plant community shape file
# December 30th 2020
# Kendra E Walters



# !!!!!!!!!     IMPORTANT: first, delete http://www.w3.org/2000/svg namespace (or whatever your file has) from the svg file  !!!!!

### Set up 

## The radius at which you want the plant composition taken for
radius <- c(0.1, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4) # set the radius for what you want plant communities defined for (1 = one meter)

## Room of ()ments
library(sp)
library(xml2)
library(tidyverse)
require(raster)
require(rgeos)
require(vegan)
require(data.table)

## Load data
file = "01_Raw_data/map_polygons__1_2_21.svg"

## Set up a few functions
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  grouped %>% group_split(.keep=F) %>% 
    rlang::set_names(rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / "))))
}

add_commas <- function(x) {
  str_split(x, " ")[[1]] %>% matrix(ncol=2, byrow=T) %>% as_tibble %>% unite(1:2, col="coord", sep=",") %>% pull(coord) %>% paste(collapse=" ") 
}


##### _______________ This be your script to make those beautiful files __________ ###########

## Importing your POLYGONS from the svg file
# note: you gotta change the script if you aren't importing polygons
pollywogs <-
  file %>% read_xml() %>% xml_find_all(".//polygon") %>% xml_attrs() %>%
  setNames(map(., ~ .x[["id"]])) %>% map(as.list) %>% map(as_tibble) %>%
  bind_rows() %>% dplyr::select(id, points) %>% mutate_at("points", function(x)
    sapply(x, add_commas)) %>%
  separate_rows(points, sep = " ") %>% separate(points, c("x", "y"), sep =
                                                  ",", convert = T) %>%
  mutate(y = -y) %>%
  named_group_split(id) %>% map(Polygon) %>% map(list) %>% map2(.y = names(.), .f =
                                                                  Polygons) %>% SpatialPolygons
row.names(pollywogs)
plot(pollywogs) # checking to see if it works

## Convert to spatial polygon dataframe
plants <- SpatialPolygonsDataFrame(pollywogs, data.frame(name=row.names(pollywogs), 
                                                         row.names=row.names(pollywogs)))
## What is a meter? 
site.coords <- as.data.frame(fread("03_Processed_data/site_locs_coordinates.tsv"))
tA1 <- site.coords[site.coords$V1 == "A1", c("x", "y")]
tA4 <- site.coords[site.coords$V1 == "A4", c("x", "y")]
tA5 <- site.coords[site.coords$V1 == "A5", c("x", "y")]
tA8 <- site.coords[site.coords$V1 == "A8", c("x", "y")]

m1 <- sqrt(sum((tA1-tA4)^2)) / 15
m2 <- sqrt(sum((tA5-tA8)^2)) / 15

tB1 <- site.coords[site.coords$V1 == "B1", c("x", "y")]
tB4 <- site.coords[site.coords$V1 == "B4", c("x", "y")]
tB5 <- site.coords[site.coords$V1 == "B5", c("x", "y")]
tB8 <- site.coords[site.coords$V1 == "B8", c("x", "y")]

m3 <- sqrt(sum((tB1-tB4)^2)) / 15
m4 <- sqrt(sum((tB5-tB8)^2)) / 15

tC1 <- site.coords[site.coords$V1 == "C1", c("x", "y")]
tC4 <- site.coords[site.coords$V1 == "C4", c("x", "y")]
tC5 <- site.coords[site.coords$V1 == "C5", c("x", "y")]
tC8 <- site.coords[site.coords$V1 == "C8", c("x", "y")]

m5 <- sqrt(sum((tC1-tC4)^2)) / 15
m6 <- sqrt(sum((tC5-tC8)^2)) / 15

one.meter <- mean(m1, m2, m3, m4, m5, m6) # 26.06
one.meter

## Prep the plot circles
xy <- site.coords[ , c("x", "y")]
names <- as.data.frame(site.coords[ , c("V1")])
sites.spt <- SpatialPointsDataFrame(coords = xy, data = names) # making it a spatial POINTS object (wild)



### Start the for loop to cycle through the different radii
for(j in 1:length(radius)) {
df.file.name <- paste0("03_Processed_data/map/grass_shrub_", radius[[j]],"_meter_composition_df.tsv")
bc.file.name <- paste0("03_Processed_data/map/grass_shrub_", radius[[j]], "_meter_bc_site_df.tsv")

## now you are ready to just run the rest of the script
site.circles <- gBuffer(sites.spt, width = radius[[j]] * one.meter, byid=TRUE)

plot(pollywogs) # adding your plant polygons for reference
plot(site.circles, add = TRUE) # checking to make sure everything looks a-okay


## INTERSECT TIME 
## to create your dataframe of plant composition by area (can be more than 100% i think) for each site
pi <- raster::intersect(plants, site.circles)
plot(plants, axes=T); plot(site.circles, add=T); plot(pi, add=T, col='red')

plant.coverage.raw <- pi@data
names(plant.coverage.raw) <- c("name", "site")

## Add in the area information
areas <- c()
for (i in 1:length(pi)) {
  areas[[i]] <- pi@polygons[[i]]@area
}

total.area <- (radius[[j]]*one.meter) * (radius[[j]]*one.meter) * 3.14159 # calculate total area
plant.coverage.raw$area.percentage <- unlist(areas) / total.area # then calculate area by percent of total

## Change those plant IDs to meaningful data
plant.key <- as.data.frame(fread("01_Raw_data/plant_composition_ID_key.txt")) # read in the key
plant.key <- subset(plant.key, select = PlantID:CroSet)
plant.coverage.raw$name <- str_extract(string = plant.coverage.raw$name, pattern = "^[a-zA-Z&&[^x]]*[0-9SMA]*") # take out the "x4" etc part of the names
plant.coverage.comm <- merge(plant.coverage.raw, plant.key, by.x = "name", by.y = "PlantID", all.x = TRUE) 
plant.coverage.comm <- plant.coverage.comm %>% mutate(across(ArtCal:CroSet, replace_na, 0)) %>% # replace all NAs with 0s
  mutate(across(ArtCal:CroSet, ~ .x * area.percentage)) # multiply the plant composition by area of intersection
plant.coverage.sum <- plant.coverage.comm %>% group_by(site) %>% summarize(across(area.percentage:CroSet, sum)) # sum for each site

## Add in rest of sites
all.site.names <- data.frame("name" = site.coords$V1) # get list of site names
plant.coverage.all <- merge(all.site.names, plant.coverage.sum, by.x = "name", by.y = "site", all.x = TRUE) # adding in sites with 100% grass

## Calculate % cover of grass for each site
plant.coverage.all <- plant.coverage.all %>% mutate(across(area.percentage:CroSet, replace_na, 0)) # turning those nasty NAs to 0s
plant.coverage.all$grass.cover <- 1 - plant.coverage.all$area.percentage # getting the percentage that should be covered in grass
plant.coverage.all$grass.cover <- ifelse(plant.coverage.all$grass.cover < 0, 0, plant.coverage.all$grass.cover) # making those over 100% to have 0% leftover space

## Adding in grassland species
plant.key <- as.data.frame(fread("01_Raw_data/plant_composition_ID_key.txt")) 
grass.vector <- subset(plant.key, subset = PlantID == "Grassland", select = ArtCal:CroSet) %>%  # get grassland data
  mutate(across(everything(), replace_na, 0)) %>% as.numeric # but no NAs!!

grass.df <- plant.coverage.all$grass.cover %o% grass.vector %>% as.data.frame # multiply grassland composition by % coverage of grass

drop <- c("area.percentage", "name", "grass.cover") # drop all non-numeric columns
row.names(plant.coverage.all) <- plant.coverage.all$name # set the row names 
row.names(grass.df) <- plant.coverage.all$name  # set the row names 
plant.cover.less <- plant.coverage.all[ , !(names(plant.coverage.all) %in% drop)] # removing those non-numeric columns

names(grass.df) <- names(plant.cover.less) # setting column names 

grass.shrub <- grass.df + plant.cover.less # adding the two matrices together == combining shrub + grass data

## making the dissimilarity  matrix
bc.grass.shrub.df <- vegdist(grass.shrub) %>% as.matrix %>% as.data.frame()

## Saving the data
write.table(bc.grass.shrub.df, file = bc.file.name, sep = "\t", row.names = TRUE)
write.table(grass.shrub, file = df.file.name, sep = "\t", row.names = TRUE)
}
#########


