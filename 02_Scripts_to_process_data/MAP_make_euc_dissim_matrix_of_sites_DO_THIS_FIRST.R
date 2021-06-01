# Make euclidean distance matrix from polylines in a SVG
# December 29th 2020
# Kendra E Walters



# !!!!!!!!!     IMPORTANT: first, delete http://www.w3.org/2000/svg namespace (or whatever your file has)    !!!!!




## Room of ()ments
require(sp)
require(xml2)
require(tidyverse)

## Load the directions
setwd("~/Google Drive/Dispersal_Ch_3/02_Data/08_Map/02_Real_Photos/10_13_20/")
file = "polygons_1_2_21.svg"

## Set up a few functions
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  grouped %>% group_split(.keep=F) %>% 
    rlang::set_names(rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / "))))
}

add_commas <- function(x) {
  str_split(x, " ")[[1]] %>% matrix(ncol=2, byrow=T) %>% as_tibble %>% unite(1:2, col="coord", sep=",") %>% pull(coord) %>% paste(collapse=" ") 
}

## Importing the weird polylines to make a euclidean distance matrix
less <- file %>% read_xml() %>% xml_find_all(".//polyline") %>% xml_attrs() %>% 
  setNames(map(., ~ .x[["id"]])) %>% map(as.list) %>% map(as_tibble) %>% 
  bind_rows() %>% dplyr::select(id, points) %>% mutate_at("points", function(x) sapply(x, add_commas)) %>% 
  separate_rows(points, sep=" ") %>% separate(points, c("x","y"), sep=",", convert=T) %>% 
  mutate(y=-y) %>% 
  named_group_split(id) %>% map(arrange, desc(y)) 

## getting site coordinates for transect A
tA <- less$TransectA
top <- data.frame("Site" = c("A1", "A2", "A3", "A4"), 
                  "x" = seq(tA[1, ]$x, tA[2, ]$x, length.out = 16)[c(1, 9, 13, 15)], 
                  "y" = seq(tA[1, ]$y, tA[2, ]$y, length.out = 16)[c(1, 9, 13, 15)])

bottom <- data.frame("Site" = c("A5", "A6", "A7", "A8"), 
                     "x" = seq(tA[2, ]$x, tA[3, ]$x, length.out = 16)[c(2, 4, 8, 16)], 
                     "y" = seq(tA[2, ]$y, tA[3, ]$y, length.out = 16)[c(2, 4, 8, 16)])
tA.sites <- rbind(top, bottom)

## getting site coordinates for transect B
tB <- less$TransectB
top <- data.frame("Site" = c("B1", "B2", "B3", "B4"), 
                  "x" = seq(tB[1, ]$x, tB[2, ]$x, length.out = 16)[c(1, 9, 13, 15)], 
                  "y" = seq(tB[1, ]$y, tB[2, ]$y, length.out = 16)[c(1, 9, 13, 15)])

bottom <- data.frame("Site" = c("B5", "B6", "B7", "B8"), 
                     "x" = seq(tB[2, ]$x, tB[3, ]$x, length.out = 16)[c(2, 4, 8, 16)], 
                     "y" = seq(tB[2, ]$y, tB[3, ]$y, length.out = 16)[c(2, 4, 8, 16)])
tB.sites <- rbind(top, bottom)

## getting site coordinates for transect C
tC <- less$TransectC
top <- data.frame("Site" = c("C1", "C2", "C3", "C4"), 
                  "x" = seq(tC[1, ]$x, tC[2, ]$x, length.out = 16)[c(1, 9, 13, 15)], 
                  "y" = seq(tC[1, ]$y, tC[2, ]$y, length.out = 16)[c(1, 9, 13, 15)])

bottom <- data.frame("Site" = c("C5", "C6", "C7", "C8"), 
                     "x" = seq(tC[2, ]$x, tC[3, ]$x, length.out = 16)[c(2, 4, 8, 16)], 
                     "y" = seq(tC[2, ]$y, tC[3, ]$y, length.out = 16)[c(2, 4, 8, 16)])
tC.sites <- rbind(top, bottom)

## merging them all together
all.sites <- rbind(rbind(tA.sites, tB.sites), tC.sites)

## Checking to make sure I didn't screw something up
all.sites$Transect <- str_extract(string = all.sites$Site, pattern = "[ABC]")
ggplot(all.sites, aes(x = x, y = y, color = Transect)) + geom_point() # looks good to me! 
all.sites$Transect <- NULL # now delete 

## Make a Euclidean DISTANCE MATRIX!!
row.names(all.sites) <- all.sites$Site
all.sites$Site <- NULL

euc.dist <- as.data.frame(as.matrix(dist(all.sites, method = "euclidean")))

## Export those bad boyz
write.table(x = euc.dist, "~/Google Drive/Dispersal_Ch_3/02_Data/site_locs_euclidean_distance_matrix.tsv", sep = "\t", row.names = TRUE)
write.table(x = all.sites, "~/Google Drive/Dispersal_Ch_3/02_Data/site_locs_coordinates.tsv", sep = "\t", row.names = TRUE)
