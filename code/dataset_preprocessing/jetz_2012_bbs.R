
library(tidyverse)
library(terra)
library(ape)


#### birds ####

## tree ##

# load data from Jetz et al. 2012, and select the first tree.
trees <- read.tree("data_raw/AllBirdsEricson1.tre")
tree <- trees[[1]] %>% ladderize()


## BBS data ##

# load datasets #

s <- list.files("data_raw/BBS/50-StopData/1997ToPresent_SurveyWide/", pattern = "csv",
                full.names = T) %>%
      map(read_csv) %>%
      map(function(x) mutate(x, StateNum = as.character(StateNum))) %>%
      bind_rows() %>%
      filter(CountryNum == 840)

sp <- readr::read_fwf("data_raw/BBS/SpeciesList.txt", skip=14) %>%
      setNames(c("Seq", "AOU", "English_Common_Name", "French_Common_Name", "Spanish_Common_Name", "ORDER", "Family", "Genus", "Species")) %>%
      select(AOU, Genus, Species)

routes <- read_csv("data_raw/BBS/routes.csv")


# merge and format #

ra <- function(a, b) round(a/b)*b
d <- s %>% select(RouteDataID:AOU) %>%
      mutate(occ = s %>% select(Stop1:Stop50) %>% rowSums()) %>%
      left_join(routes) %>%
      group_by(AOU, Longitude, Latitude) %>%
      summarize(occ = sum(occ)) %>%
      filter(occ > 0) %>%
      ungroup() %>%
      left_join(sp) %>%
      mutate(binomial = paste(Genus, Species)) %>%
      mutate(x = ra(Longitude, 2),
             y = ra(Latitude, 2)) %>%
      filter(y <= 50,
             !str_detect(binomial, "sp\\."),
             !str_detect(binomial, " / "),
             !str_detect(binomial, " x "))
g <- d %>%
      group_by(x, y) %>%
      summarize(z = length(unique(binomial))) %>%
      rast() %>%
      setValues(NA)
r <- d %>%
      split(.$binomial) %>%
      map(function(f) rasterize(as.matrix(select(f, x, y)), g)) %>%
      rast()


# harmonize and export #

tree$tip.label <- str_replace(tree$tip.label, "_", " ")
species <- intersect(names(r), tree$tip.label)

r <- r[[species]]
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

writeRaster(r, "data/jetz_2012/bbs_occ.tif")
write.nexus(tree, file = "data/jetz_2012/jetz_2012_tree0001.nex")
