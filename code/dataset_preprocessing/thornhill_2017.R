

library(tidyverse)
library(ape)


# function to resolve polytomies by matching a reference tree topology
multi2di_target <- function(x, y){
      require(phytools)
      require(ape)

      while(! is.binary(x)){

            # identify the nodes with polytomies
            node <- x$edge %>% as.data.frame() %>% group_by(V1) %>% count() %>% filter(n > 2)
            node <- node$V1[1]

            x <- splitTree(x, list(node=node,
                                   bp=x$edge.length[which(x$edge[,2] == node)]))
            names(x) <- c("backbone", "clade")

            # resolve polytomies to match reference clade (inefficient trial & error method)
            target <- extract.clade(y, getMRCA(y, x$clade$tip.label))
            graft <- multi2di(x$clade, random=T)
            while(! all.equal(graft, target, use.edge.length=F)) graft <- multi2di(x$clade, random=T)

            # reattach modified clade to backbone tree
            x <- bind.tree(x$backbone, graft, which(x$backbone$tip.label == "NA"))

      }

      return(x)
}


# load phylogeny
chrono <- read.nexus("data_raw/thornhill_2017/doi_10.6078_D1VD4P__v3/California_clades_fully_dated.nex")
phylo <- read.nexus("data_raw/thornhill_2017/doi_10.6078_D1VD4P__v3/Californian_clades_tree_final.nex")

# the chronogram has polytomies that the phylogram does not.
# fix this, and then rotate so that edge matrices correspond.
chrono <- multi2di_target(chrono, phylo)
chrono <- rotateConstr(chrono, phylo$tip.label)
phylo <- rotateConstr(phylo, phylo$tip.label)

# double check that all edges correspond
for(i in 1:nrow(chrono$edge)){
      if(chrono$edge[i,2] > length(chrono$tip.label)){
            a <- extract.clade(chrono, chrono$edge[i,2])$tip.label
            b <- extract.clade(phylo, phylo$edge[i,2])$tip.label
            if(length(c(setdiff(a, b), setdiff(b, a))) > 0) stop("problem")
      }else{
            a <- chrono$tip.label[chrono$edge[i,2]]
            b <- phylo$tip.label[phylo$edge[i,2]]
            if(a != b) stop("problem")
      }
}

# normalize edge lengths
chrono$edge.length <- chrono$edge.length / sum(chrono$edge.length)
phylo$edge.length <- phylo$edge.length / sum(phylo$edge.length)

# load occurrence points
occ <- read_csv("data_raw/thornhill_2017/California_Clades_clean_All_final.csv",
                col_select = c(x = longitude, y = latitude, otu = clade))

# clean taxon names
chrono$tip.label <- str_replace_all(chrono$tip.label, "\\+", "__")
phylo$tip.label <- str_replace_all(phylo$tip.label, "\\+", "__")
occ$otu <- str_replace_all(occ$otu, "\\+", "__")

# select overlapping species
otu <- intersect(chrono$tip.label, unique(occ$otu))
occ <- occ[occ$otu %in% otu,]

# load and aggregate data on current land protection status
# (from Kling et al. 2018)
library(raster)
rnd <- function(x, y) round(x/y)*y
grid_res <- 0.5
p <- raster("data_raw/kling_2018_protection_status.tif") %>%
      projectRaster(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>%
      rasterToPoints() %>% as.data.frame() %>% as_tibble() %>%
      mutate(x = rnd(x, grid_res), y = rnd(y, grid_res)) %>%
      group_by(x, y) %>%
      summarize(p = mean(kling_2018_protection_status)) %>%
      rasterFromXYZ()


# export
write_csv(occ, "data/thornhill_2017/thornhill_2017_occ.tif")
write.nexus(chrono, file = "data/thornhill_2017/thornhill_2017_chronogram.nex")
write.nexus(phylo, file = "data/thornhill_2017/thornhill_2017_phylogram.nex")
writeRaster(p, "data/thornhill_2017/kling_2018_protection_status.tif", overwrite = T)
