# Print some descriptive stats about a graph

descriptive_stats <- function(g) {
    x <- list("order" = igraph::gorder(g),
              "size" = igraph::gsize(g))
    is_directed <- igraph::is_directed(g)
    if (is_directed(g)) {
        x[["in-degree"]] <- round(mean(igraph::degree(graph = g, mode = "in")), 3)
        x[["out-degree"]] <- round(mean(igraph::degree(graph = g, mode = "out")), 3)
        x[["mean-distance"]] <- round(igraph::mean_distance(graph = g,
                                                      directed = TRUE,
                                                      unconnected = TRUE), 3)
        x[["diameter"]] <- igraph::diameter(graph = g,
                                            directed = FALSE,
                                            unconnected = TRUE)
    } else {
        x[["degree"]] <- round(mean(igraph::degree(graph = g, mode = "all")), 3)
        x[["mean-distance"]] <- round(igraph::mean_distance(graph = g,
                                                      directed = FALSE,
                                                      unconnected = TRUE), 3)
        x[["diameter"]] <- igraph::diameter(graph = g,
                                            directed = TRUE,
                                            unconnected = TRUE)
    }
    cat(paste("Order:\t\t\t", x[["order"]], "\n"))
    cat(paste("Size:\t\t\t", x[["size"]], "\n"))
    cat(paste("Diameter:\t\t", x[["diameter"]], "\n"))
    cat(paste("Average Path Length:\t", x[["mean-distance"]], "\n"))
    if (is_directed) {
        cat(paste("Average in-degree:\t", x[["in-degree"]], "\n"))
        cat(paste("Average out-degree:\t", x[["out-degree"]], "\n"))
    } else {
        cat(paste("Average degree:\t\t", x[["degree"]], "\n"))
    }

}