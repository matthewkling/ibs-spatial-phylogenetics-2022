
library(sf)
library(tidyverse)
library(terra)
library(raster)
library(ape)

# mammal phylogeny
# from Upham et al. 2019: http://vertlife.org/data/mammals/
# (the paper has a large set of posterior trees; we'll use just one)
tree <- read.tree("data_raw/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
write.tree(tree, "data/upham_2019/upham_2019_tree0000.tre")

# terestrial mammal range polygons
# from IUCN: https://www.iucnredlist.org/resources/spatial-data-download
m <- st_read("data_raw/MAMMALS_TERRESTRIAL_ONLY/")

# grid template
gr <- expand_grid(x = -179:180 - .5,
                  y = -89:90 - .5,
                  z = 1) %>%
      rast() %>%
      raster()
g <- gr %>%
      rasterToPolygons() %>%
      st_as_sf()
st_crs(g) <- st_crs(m)

# gridded presences
sf_use_s2(FALSE) # switching off spherical geometry to prevent edge-crossing error
x <- st_intersects(g, m, sparse = F) # faster to use sf than raster::rasterize

# format output
xd <- x %>%
      t() %>%
      as.data.frame() %>%
      setNames(paste0("cell_", 1:nrow(x))) %>%
      mutate(species = m$binomial) %>%
      group_by(species) %>%
      summarize_all(any) %>%
      ungroup() %>%
      column_to_rownames("species") %>%
      as.matrix() %>%
      t()
xr <- apply(xd, 2, function(z) setValues(gr, z)) %>%
      map(rast) %>%
      rast() %>%
      setNames(colnames(xd)) %>%
      aggregate(3, fun = max)
writeRaster(xr, "data/upham_2019/iucn_mammals.tif", overwrite = T)

# climate data
clim <- list.files("~/climate_data/CHELSA/v2/", full.names = T) %>%
      rast() %>%
      setNames(c("tmean", "precip", "tmax", "tmin")) %>%
      aggregate(400) %>%
      resample(xx[[1]])
writeRaster(clim, "data/upham_2019/chelsa_climate.tif", overwrite = T)


