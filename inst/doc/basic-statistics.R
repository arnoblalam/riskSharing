## ----short description of data-------------------------------------------
library(riskSharing)
data("nyakatoke")
knitr::kable(head(nyakatoke[,1:4]))

## ----households pairs appear twice---------------------------------------
knitr::kable(subset(nyakatoke[,1:4], (hh1 == 1 & hh2 == 76) | (hh1 == 76 & hh2 == 1)))

## ----head all columns----------------------------------------------------
knitr::kable(head(nyakatoke))

## ----create networks-----------------------------------------------------
# For the desire to link data, we create an "edgelist" that igraph can read.
# This is a list of households where the first household in the pair mentioned
# someone from household 2 in their survey response.
edgelist <- as.matrix(nyakatoke[nyakatoke$willingness_link1 == 1, 1:2])
g.directed <- igraph::graph_from_data_frame(edgelist)

# For the underreporting model, we create a graph from the dataframe
# First we keep only the rows where at least one houshold said the other is in
# their network.
underreporting.df <- 
    nyakatoke[nyakatoke$willingness_link1 == 1 | nyakatoke$willingness_link2 == 1,]
# We then convert this into an undirected graph
g.underreporting <- igraph::graph_from_data_frame(underreporting.df, 
                                                  directed = FALSE)
# We then remove any multiple edges
g.underreporting <- igraph::simplify(g.underreporting)

# For the overreporting model we only keep those rows where both households
# listed the other in their network
overreporting.df <- 
    nyakatoke[nyakatoke$willingness_link1 == 1 & nyakatoke$willingness_link2 == 1,]
# We then create an graph
g.overreporting <- igraph::graph_from_data_frame(overreporting.df, directed = FALSE)
# From examining the list of vertices, we note that households 7, 30 32, 36, 44,
# 65, 84, 88, 91, 96, 107, 110, 116, 117, 118 and 119 don't show up in the 
# graph. This is because there is no case where these households list another 
# household in their support network and they do the same.  We dd back these 
# vertices as isolated vertices
missing.vertices <- c(7, 30, 32, 46, 44, 65, 84, 88, 91, 96, 107, 110, 
                      116, 117, 118, 119)
g.overreporting <- (Reduce(f = function(x, y) {y + igraph::vertex(x)},
       x = missing.vertices,
       init = g.overreporting,
       right = TRUE))
# We remove multiple edges
g.overreporting <- igraph::simplify(g.overreporting)

## ----descriptive stats directed------------------------------------------
descriptive_stats(g.directed)

## ----descriptive stats under---------------------------------------------
descriptive_stats(g.underreporting)

## ----descriprive stats over----------------------------------------------
descriptive_stats(g.overreporting)

