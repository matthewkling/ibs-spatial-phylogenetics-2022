
# function to clean a spatial phylogenetic dataset
clean_dataset <- function(x){

      # remove spaces from taxon names
      x$tree$tip.label <- str_replace(x$tree$tip.label, " ", "_")
      colnames(x$comm) <- str_replace(colnames(x$comm), " ", "_")

      # convert non-finite values to formal absences
      x$comm[!is.finite(x$comm)] <- 0

      # remove taxa with no occurrences
      ghost <- colSums(x$comm, na.rm = T) == 0
      if(any(ghost)){
            message("removing ", sum(ghost), " taxa with no occurrences")
            x$comm <- x$comm[, !ghost]
      }

      # remove non-overlapping taxa
      taxa <- intersect(x$tree$tip.label, colnames(x$comm)[colSums(x$comm, na.rm = T) > 0])
      if(length(taxa) == 0) stop("no overlapping taxon names")
      prune <- setdiff(x$tree$tip.label, taxa)
      if(length(prune > 0)){
            message("removing ", length(prune), " tips from tree")
            x$tree <- drop.tip(x$tree, prune)
      }
      drop <- setdiff(colnames(x$comm), taxa)
      if(length(drop > 0)){
            message("removing ", length(drop), " taxa from community matrix")
            x$comm <- x$comm[, taxa]
      }

      # remove sites with no taxa
      occupied <- rowSums(x$comm) > 0
      if(sum(occupied > 0)) message("removing ", sum(occupied), " sites with no occurrences")
      x$comm <- x$comm[occupied,]
      x$xy <- x$xy[occupied,]

      # name rows
      rownames(x$comm) <- paste0("cell_", 1:nrow(x$comm))

      return(x)
}

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
