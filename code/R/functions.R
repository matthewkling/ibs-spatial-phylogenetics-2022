
parentProb <- function(x) 1 - prod(1 - x)

build_clade_range <- function(e, phylo, sxt){
      node <- phylo$edge[e,2]
      if(node <= length(phylo$tip.label)){
            otu <- phylo$tip.label[node]
            prob <- sxt[,otu]
      } else{
            clade <- extract.clade(phylo, node)
            otu <- clade$tip.label
            prob <- apply(sxt[,otu], 1, parentProb)
      }
      return(prob)
}

build_clade_ranges <- function(tree, tip_occs){
      sapply(1:nrow(tree$edge), build_clade_range, phylo=tree, sxt=tip_occs)
}

plot_edge_connect <- function(tree, tips, hl = "orange", bg = "black", ...){
      clr <- rep(bg, length(tree$edge.length))
      clr[which.edge(tree, tips)] <- "orange"
      plot(tree, edge.color = clr, ...)
}

plot_edge_numeric <- function(tree, var, ...){
      pal <- colorRampPalette(c("black", "blue", "red", "orange"))(20)
      clr <- pal[cut(branch_length, 20)]
      if(sd(var[is.finite(var)]) == 0) clr <- "black"
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
