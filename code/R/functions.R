
# a set of functions to construct a community matrix
# with a row for every branch of a phylogeny
parentProb <- function(x) 1 - prod(1 - x)
build_clade_range <- function(e, phylo, comm){
      node <- phylo$edge[e,2]
      if(node <= length(phylo$tip.label)){
            otu <- phylo$tip.label[node]
            prob <- comm[,otu]
      } else{
            clade <- extract.clade(phylo, node)
            otu <- clade$tip.label
            prob <- apply(comm[,otu], 1, parentProb)
      }
      return(prob)
}
build_clade_ranges <- function(tree, comm){
      sapply(1:nrow(tree$edge),
             build_clade_range, phylo=tree, comm=comm)
}

# helpher function to plot a phylogeny
phyplot <- function(tree,
                   connect = NULL, hl = "orange", bg = "black",
                   value = NULL, col = c("black", "blue", "red", "orange"),
                   ...){

      clr <- rep(bg, length(tree$edge.length))
      if(!is.null(connect)){
            clr[which.edge(tree, connect)] <- hl
      }

      if(!is.null(value)){
            n <- min(c(20, length(tree$edge.length)))
            pal <- colorRampPalette(col)(n)
            clr <- pal[cut(value, n)]
            if(sd(value[is.finite(value)]) == 0) clr <- "black"
      }

      plot(tree, edge.color = clr, ...)
}

# helper function to plot a map
carto <- function(d, v){
      # md <- map_data("state", "california")
      # md <- map_data("state")
      p <- ggplot() +
            geom_tile(data = d,
                      aes_string("x", "y", fill = v)) +
            # geom_path(data = md,
            #           aes(long, lat, group = group),
            #           color = "white", alpha = .5) +
            xlim(range(d$x) + diff(range(d$x))/20 * c(-1, 1)) +
            ylim(range(d$y) + diff(range(d$y))/20 * c(-1, 1)) +
            coord_fixed(ratio = 1.2) +
            theme_void() +
            theme(legend.position = "top")
      # theme(legend.position = c(.8, .95),
      #       legend.justification = c(1, 1))
      if(inherits(d[[v]], "numeric")) p <- p + scale_fill_viridis_c()
      p
}
