
library(tidyverse)
library(ape)

# load points and phylogeny from mishler et al. 2014
d <- read_csv("data_raw/Point_distribution_Australian_Phylogenetic_Diversity_Acacia.csv")
p <- read.nexus("data_raw/1_1363828941_Acacia.nex")

# format tip labels
p$tip.label <- p$tip.label %>%
      tolower() %>%
      str_remove("acacia") %>%
      str_replace_all("_|\\.", " ") %>%
      str_trim() %>%
      str_split(" ", n = 2, simplify = T) %>%
      "["( , 1)

# intersect taxa and export #
otu <- intersect(p$tip.label, unique(d$taxa))
d %>%
      filter(taxa %in% otu) %>%
      write_csv("data/mishler_2014/mishler_2014_acacia_points.csv")
p %>%
      drop.tip(setdiff(p$tip.label, otu)) %>%
      write.tree("data/mishler_2014/mishler_2014_acacia.tre")
