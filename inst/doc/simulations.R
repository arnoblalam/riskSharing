## ------------------------------------------------------------------------
library(riskSharing)
data("nyakatoke")

edgelist <- as.matrix(nyakatoke[nyakatoke$willingness_link1 == 1, 1:2])
g.directed <- igraph::graph_from_data_frame(edgelist)

underreporting.df <- 
    nyakatoke[nyakatoke$willingness_link1 == 1 | nyakatoke$willingness_link2 == 1,]
g.underreporting <- igraph::graph_from_data_frame(underreporting.df, 
                                                  directed = FALSE)
g.underreporting <- igraph::simplify(g.underreporting)

overreporting.df <- 
    nyakatoke[nyakatoke$willingness_link1 == 1 & nyakatoke$willingness_link2 == 1,]
g.overreporting <- igraph::graph_from_data_frame(overreporting.df, 
                                                 directed = FALSE)
missing.vertices <- c(7, 30, 32, 46, 44, 65, 84, 88, 91, 96, 107, 110, 
                      116, 117, 118, 119)
g.overreporting <- (Reduce(f = function(x, y) {y + igraph::vertex(x)},
       x = missing.vertices,
       init = g.overreporting,
       right = TRUE))

## ----simulation parameters-----------------------------------------------
sim.size <- 1000
sim.order <- 119
size.directed <- igraph::gsize(g.directed)
size.underreporting <- igraph::gsize(g.underreporting)
size.overreporting <- igraph::gsize(g.overreporting)
avg.degree.directed <- mean(igraph::degree(g.directed, mode = "out"))
avg.degree.underreporting <- mean(igraph::degree(g.underreporting))
avg.degree.overreporting <- mean(igraph::degree(g.overreporting))

## ----E-R-M underreporting simulation-------------------------------------
underreporting.sim <- replicate(n = sim.size,
                                igraph::sample_gnm(sim.order, size.underreporting), 
                                simplify = FALSE)

## ----underreporting size order and degree--------------------------------
# The order of all the graphs should be sim.order = 119 and the size should be 
# size.underreporting = 490
all(sapply(underreporting.sim, igraph::gorder) == sim.order)
all(sapply(underreporting.sim, igraph::gsize) == size.underreporting)
all(sapply(underreporting.sim, function (x) mean(igraph::degree(x)))
        ==  2*size.underreporting/sim.order)

## ----underreporting clustering coefficient, fig.width=6, fig.height=5----
underreporting.cc <- sapply(underreporting.sim, 
                            function(x) igraph::transitivity(x))
summary(underreporting.cc)
round(quantile(underreporting.cc, c(0.025, 0.975)), 3)
h <- hist(underreporting.cc, plot = FALSE)
h$density <- h$counts/sum(h$counts)
plot(h, 
     freq = FALSE,
     xlab = "Clustering Coefficient",
     main = "Distribtuion of Global Clustering Coefficients
     in Erdos-Renyi simulations")

## ----underreporting local clustering, fig.width=6, fig.height=5----------
underreporting.local.cc <- sapply(underreporting.sim, function(g) {
    mean(igraph::transitivity(graph = g, type = "local", isolates = "zero"))
})

round(quantile(underreporting.local.cc, c(0.025, 0.975)), 3)
h <- hist(underreporting.local.cc, plot = FALSE)
h$density <- h$counts/sum(h$counts)
plot(h, 
     freq = FALSE,
     xlab = "Local Clustering Coefficient",
     main = "Distribtuion of Local Clustering Coefficients
     in Erdos-Renyi simulations")

## ----overreporting simulation--------------------------------------------
overreporting.sim <- replicate(n = sim.size,
                                igraph::sample_gnm(sim.order, size.overreporting), 
                                simplify = FALSE)

