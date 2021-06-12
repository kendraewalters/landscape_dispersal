# Script to create plant composition data for the grassland and shrubland ecosystems
# June 9th 2021
# Kendra E Walters

# Room of ()ments
require(data.table)
library(sp)
library(xml2)
library(tidyverse)
require(raster)
require(rgeos)
require(vegan)


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


## Setting up what is a meter
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


## Importing your POLYGONS from the svg file
# note: you gotta change the script if you aren't importing polygons
plants <-
  file %>% read_xml() %>% xml_find_all(".//polygon") %>% xml_attrs() %>%
  setNames(map(., ~ .x[["id"]])) %>% map(as.list) %>% map(as_tibble) %>%
  bind_rows() %>% dplyr::select(id, points) %>% mutate_at("points", function(x)
    sapply(x, add_commas)) %>%
  separate_rows(points, sep = " ") %>% separate(points, c("x", "y"), sep =
                                                  ",", convert = T) %>%
  mutate(y = -y) %>%
  named_group_split(id) %>% 
  map(Polygon) %>% 
  map(list) %>% 
  map2(.y = names(.), .f = Polygons) %>% 
  SpatialPolygons %>% 
  SpatialPolygonsDataFrame(data.frame(name=row.names(.), 
                                      row.names=row.names(.)))


## Create the lat lons for our two rectangles (sort of), one for grassland, one for shrubland
# Making the boundary points between the two ecosystems, including the 4m max radius
boundary_points <- fread("03_Processed_data/site_locs_coordinates.tsv", data.table = FALSE) %>% 
  filter(grepl("[ABC][45]", V1)) %>% 
  mutate(Transect = str_extract(V1, "[ABC]")) %>% 
  group_by(Transect) %>% 
  summarize(across(c(x,y), .fns = mean)) %>% 
  mutate(x_plus_4m = x + 4*one.meter, 
         x_minus_4m = x - 4*one.meter) %>% 
  pivot_longer(contains("x"), names_to = "x_position", values_to = "x_coord") %>% 
  filter(Transect == "A" & x_position == "x_plus_4m" | 
           Transect == "B" & x_position == "x" |
           Transect == "C" & x_position == "x_minus_4m") %>% 
  select(!c(x_position, Transect)) %>% 
  rename("x" = x_coord) %>% 
  mutate(n = 2) %>% 
  uncount(n) %>% 
  mutate(ecosystem = ifelse(duplicated(.) == TRUE, "grassland", "shrubland")) # just trying to assign same values to shrubland and grassland

# Making the boundary points of our transects, including the 4m max radius
# And then merging with previous dataset and creating the polygon rectangles
ecosystem.polygon <- fread("03_Processed_data/site_locs_coordinates.tsv", data.table = FALSE) %>% 
  filter(grepl("[ABC][18]", V1)) %>% 
  mutate(x_plus_4m = x + 4*one.meter, 
         x_minus_4m = x - 4*one.meter, 
         y_plus_4m = y + 4*one.meter, 
         y_minus_4m = y - 4*one.meter) %>% 
  pivot_longer(contains("x"), names_to = "x_position", values_to = "x_coord") %>% 
  pivot_longer(contains("y"), names_to = "y_position", values_to = "y_coord") %>% 
  filter(V1 == "A1" & x_position == "x_plus_4m" & y_position == "y_plus_4m" | 
           V1 == "A8" & x_position == "x_plus_4m" & y_position == "y_minus_4m" |
           V1 == "B1" & x_position == "x" & y_position == "y_plus_4m" |
           V1 == "B8" & x_position == "x" & y_position == "y_minus_4m" |
           V1 == "C1" & x_position == "x_minus_4m" & y_position == "y_plus_4m" |
           V1 == "C8" & x_position == "x_minus_4m" & y_position == "y_minus_4m") %>% 
  mutate(ecosystem = ifelse(y_position == "y_plus_4m", "shrubland", "grassland")) %>% 
  select(!c(x_position, y_position, V1)) %>% 
  rename("x" = x_coord, 
         "y" = y_coord) %>% 
  full_join(boundary_points) %>% 
  mutate(polygon_order = 1:nrow(.), 
         polygon_order = ifelse(polygon_order == 7, 11,
                                ifelse(floor(x) == 428 & ecosystem == "shrubland", 7, 
                                       ifelse(polygon_order == 8, 12, 
                                              ifelse(floor(x) == 428 & ecosystem == "grassland", 8, polygon_order))))) %>% 
  arrange(polygon_order) %>% # this was a stupid way to just order them correctly so the polygon lines wouldn't cross (probably a better way to do this??)
  select(!polygon_order) %>% 
  named_group_split(ecosystem) %>% 
  map(Polygon) %>% 
  map(list) %>%
  map2(.y = names(.), .f = Polygons) %>%
  SpatialPolygons() %>% 
  SpatialPolygonsDataFrame(data.frame(name=row.names(.), 
                                      row.names=row.names(.)))

plot(plants)
plot(ecosystem.polygon, add = TRUE) # checking to see if it works


## Intersecting our two polygons
pi <- raster::intersect(plants, ecosystem.polygon)
plot(plants, axes=T); plot(ecosystem.polygon, add=T); plot(pi, add=T, col='red') # checking to see if it worked

areas <- c()
for (i in 1:length(pi)) {
  areas[[i]] <- pi@polygons[[i]]@area
}

plant.key <- as.data.frame(fread("01_Raw_data/plant_composition_ID_key.txt")) %>% 
  select(PlantID:CroSet)

plant.coverage.ecosystem <- pi@data %>% 
  rename("name" = name.1, 
         "ecosystem" = name.2) %>% 
  mutate(area = unlist(areas)) %>% 
  add_row(name = "Grassland", ecosystem = "shrubland", # adding new row for grassland (which I didn't code directly but left as total area - everything else)
          area = ecosystem.polygon@polygons$shrubland@area - # total area of shrubland polygon - area already filled with plants or bare ground
            aggregate(.$area, list(.$ecosystem), sum)$x[2]) %>% 
  add_row(name = "Grassland", ecosystem = "grassland", 
          area = ecosystem.polygon@polygons$grassland@area - # total area of grassland polygon - area already filled with plants or bare ground
            aggregate(.$area, list(.$ecosystem), sum)$x[1]) %>% 
  left_join(plant.key, by = c("name" = "PlantID")) %>% 
  mutate(area_m2 = area / (one.meter * one.meter)) %>% 
  mutate(across(ArtCal:CroSet, replace_na, 0)) %>% # replace all NAs with 0s
  mutate(across(ArtCal:CroSet, ~ .x * area_m2)) %>% # multiply the plant composition by area of intersection
  group_by(ecosystem) %>% summarize(across(ArtCal:area_m2, sum)) %>% # sum for each ecosystem
  mutate(BareGround = area_m2 - rowSums(across(ArtCal:CroSet)))


## Export table
write_tsv(plant.coverage.ecosystem, "03_Processed_data/plant_composition_m2_in_each_ecosystem_site.tsv")







