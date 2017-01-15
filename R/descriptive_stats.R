# Print some descriptive stats about a graph

descriptive_stats <- function(g) {
    x <- list("order" = igraph::gorder(g),
              "size" = igraph::gsize(g))
    is_directed <- igraph::is_directed(g)
    if (is_directed(g)) {
        x[["in-degree"]] <- mean(igraph::degree(graph = g, mode = "in"))
        x[["out-degree"]] <- mean(igraph::degree(graph = g, mode = "out"))
        x[["mean-distance"]] <- igraph::mean_distance(graph = g,
                                                      directed = TRUE,
                                                      unconnected = TRUE)
        x[["diameter"]] <- igraph::diameter(graph = g,
                                            directed = FALSE,
                                            unconnected = TRUE)
    } else {
        x[["degree"]] <- mean(igraph::degree(graph = g, mode = "all"))
        x[["mean-distance"]] <- igraph::mean_distance(graph = g,
                                                      directed = FALSE,
                                                      unconnected = TRUE)
        x[["diameter"]] <= igraph::diameter(graph = g,
                                            directed = TRUE,
                                            unconnected = TRUE)
    }
}