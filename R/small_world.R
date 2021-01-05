#### Project title: The structure of risk-sharing networks
#### Authors: H. Henderson and A. Alam
#### Last updated: 1/1/2021
#### Purpose: Creates results for Section 3.1 of the paper

# Set up and read in data
setwd("~/Documents/Articles/Risk-Sharing-Networks/Code")
nyakatoke <- read.csv("Nyakatoke.csv")
library(igraph)
library(doParallel)
library(parallel)
library(foreach)

# Register parallel backends
cores <- detectCores()
registerDoParallel(cores = cores)

# Function to calculate global clustering coefficient for directed and undirected networks
clustering <- function(network){
  n <- igraph::vcount(network)
  num <- 0       
  den <- 0
  for (i in 1:n){       # Loop through each node i, finding i's neighbors (excluding i)
    n1 <- unlist(igraph::ego(network, order=1, nodes=i, mode="out"))
    n1 <- n1[n1 != i]
    for (j in n1){      # Loop through i's neighbors j, finding j's neighbors (excluding i and j)
      n2 <- unlist(igraph::ego(network, order=1, nodes=j, mode="out"))
      n2 <- n2[n2 != i & n2 != j]
      den <- den + length(n2)   # Add number of j's neighbors to denominator
      for (k in n2){    # Loop through each of j's neighbors k
        n3 <- unlist(igraph::ego(network, order=1, nodes=k, mode="out"))  
        num <- num + is.element(i, n3)  # Increment numerator if i is a neighbor of k
      }
    }
  }
  cc <- num/den
  return(cc)
}

# Desire-to-link network
edgelist <- as.matrix(nyakatoke[nyakatoke$willingness_link1 == 1, 1:2])
g.desire <- igraph::graph_from_data_frame(edgelist)

# Underreporting network
underreporting.df <- 
  nyakatoke[nyakatoke$willingness_link1 == 1 | nyakatoke$willingness_link2 == 1,]

# Read underreporting network into igraph and remove redundant links
g.underreporting <- igraph::graph_from_data_frame(underreporting.df, directed = FALSE)
g.underreporting <- igraph::simplify(g.underreporting)

# Overreporting network
overreporting.df <- 
  nyakatoke[nyakatoke$willingness_link1 == 1 & nyakatoke$willingness_link2 == 1,]

# Read overreporting network into igraph
g.overreporting <- igraph::graph_from_data_frame(overreporting.df, directed = FALSE)

# From examining the list of vertices, we find that households 7, 30 32, 36, 44,
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

# Remove multiple edges from overreporting network
g.overreporting <- igraph::simplify(g.overreporting)

# Select network to use: g.underreporting, g.overreporting, or g.desire
network <- g.desire
if (identical_graphs(network, g.underreporting)){name <- "Underreporting"}
if (identical_graphs(network, g.overreporting)){name <- "Overreporting"}
if (identical_graphs(network, g.desire)){name <- "Desire-to-link"}
directed <- igraph::is.directed(network)

# Number of nodes and links
n <- igraph::vcount(network)
links <- igraph::ecount(network)

# Average degree
deg <- round(mean(igraph::degree(network, mode="in")), 2)

# Percent in giant component
comp <- max(igraph::components(network)$csize)
comp <- round(100*(comp/n), 2)

# Average distance
dist <- round(igraph::mean_distance(network), 2)

# Clustering coefficient
# Note: For undirected networks, our clustering function and igraph's (global) transitivity
# function give the same result. For those networks, we use igraph's function because it 
# is faster. For directed networks, our clustering function is different from igraph's and
# we use our function because it is appropriate for directed networks.
clus <- ifelse(name=="Desire-to-link",
  round(clustering(network), 2),  
  round(igraph::transitivity(network, type="global"), 2))

# Simulate corresponding random network
reps <- 10000
deg.sim <- numeric(reps)
comp.sim <- numeric(reps)
dist.sim <- numeric(reps)
clus.sim <- numeric(reps)
set.seed(12345)
results <- foreach(i = 1:reps) %dopar% {
  ntw <- igraph::sample_gnm(n, links, directed = directed)
  deg.sim <- mean(igraph::degree(ntw, mode="in"))
  comp.sim <- max(igraph::components(ntw)$csize)/n
  dist.sim <- igraph::mean_distance(ntw)
  clus.sim <- ifelse(name=="Desire-to-link",
          clustering(ntw),  
          igraph::transitivity(ntw, type="global")
  )
  res <- c(deg.sim, comp.sim, dist.sim, clus.sim)
}
results <- as.data.frame(do.call(rbind, results)) # Aggregate results in data frame
deg.sim <- results[, 1]  # Unpack results
comp.sim <- results[, 2]
dist.sim <- results[, 3]
clus.sim <- results[, 4]
deg.avg <- round(mean(deg.sim), 2)   # Calculate averages
comp.avg <- round(100*mean(comp.sim), 2)
dist.avg <- round(mean(dist.sim), 2)
clus.avg <- round(mean(clus.sim), 2)
comp.sd <- round(sd(comp.sim), 2) # Calculate SDs
dist.sd <- round(sd(dist.sim), 2)
clus.sd <- round(sd(clus.sim), 2)

# Print results
cat(
  paste0("Network: ", name),
  paste0("Number of nodes: ", n),
  paste0("Number of links: ", links),
  paste0("Average degree: ", deg),
  paste0("Percent in giant component: ", comp),
  paste0("Average distance: ", dist),
  paste0("Clustering coefficient: ", clus),
  "",
  paste0("Network: Random"),
  paste0("Number of nodes: ", n),
  paste0("Number of links: ", links),
  paste0("Average degree: ", deg.avg),
  paste0("Percent in giant component: ", comp.avg, " (", comp.sd, ")"),
  paste0("Average distance: ", dist.avg, " (", dist.sd, ")"),
  paste0("Clustering coefficient: ", clus.avg, " (", clus.sd, ")"),
  sep="\n")

# stopCluster(cl)



