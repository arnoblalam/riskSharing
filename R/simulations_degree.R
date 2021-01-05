#### Project title: The structure of risk-sharing networks
#### Authors: H. Henderson and A. Alam
#### Last updated: 1/1/2021
#### Purpose: Creates results for Section 4.1 of the paper

# Set up and read in data
setwd("~/Documents/Articles/Risk-Sharing-Networks/Code")
library(igraph)
source("sim_functions_degree.R")
nyakatoke <- read.csv("Nyakatoke.csv")
set.seed(12345)

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

# Define x-axis
x <- (0:119/119)

### Simulate random networks and calculate average distance
reps <- 10000
under.dist <- numeric(reps)
over.dist <- numeric(reps)
desire.dist <- numeric(reps)
for (i in 1:reps){
  under.rand <- jr(119, 490, pr=1, pn=0, directed=FALSE)
  over.rand <- jr(119, 140, pr=1, pn=0, directed=FALSE)
  desire.rand <- jr(119, 630, pr=1, pn=0, directed=TRUE)
  
  res.under <- mean_distance(under.rand, directed=FALSE)
  res.over <- mean_distance(over.rand, directed=FALSE)
  res.desire <- mean_distance(desire.rand, directed=TRUE)
  
  under.dist[i] <- res.under
  over.dist[i] <- res.over
  desire.dist[i] <- res.desire
}
print(mean(under.dist))
print(mean(over.dist))
print(mean(desire.dist))

### Simulations for underreporting network

# Robustness to random attack
under.random <- degree.sims(g.underreporting, "random", 10000) 
under.random <- rbind(c(1, 1, 1), under.random, c(0, 0, 0), c(0, 0, 0)) 
under.random <- cbind(x, under.random)

# Robustness to targeted attack
under.targeted <- degree.sims(g.underreporting, "targeted", 10000)
under.targeted <- rbind(c(1, 1, 1), under.targeted, c(0, 0, 0), c(0, 0, 0)) 
under.targeted <- cbind(x, under.targeted)

### Simulations for overreporting network

# Robustness to random attack
over.random <- degree.sims(g.overreporting, "random", 10000) 
over.random <- rbind(c(1, 1, 1), over.random, c(0, 0, 0), c(0, 0, 0)) 
over.random <- cbind(x, over.random)

# Robustness to targeted attack
over.targeted <- degree.sims(g.overreporting, "targeted", 10000)
over.targeted <- rbind(c(1, 1, 1), over.targeted, c(0, 0, 0), c(0, 0, 0)) 
over.targeted <- cbind(x, over.targeted)

### Simulations for desire-to-link network

# Robustness to random attack
desire.random <- degree.sims(g.desire, "random", 10000) 
desire.random <- rbind(c(1, 1, 1), desire.random, c(0, 0, 0), c(0, 0, 0)) 
desire.random <- cbind(x, desire.random)

# Robustness to targeted attack
desire.targeted <- degree.sims(g.desire, "targeted", 10000)
desire.targeted <- rbind(c(1, 1, 1), desire.targeted, c(0, 0, 0), c(0, 0, 0)) 
desire.targeted <- cbind(x, desire.targeted)

### Graph for robustness to random attack
m <- matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4), nrow=10, ncol=1)
layout(m)
par(mar = c(4, 4, 3, 1))

plot(c(0, 1), c(0, 1), type="n", xlab="Proportion of nodes removed", # Plot results
     ylab="Relative size of giant component", main="(a) Underreporting")
lines(x, under.random$results.rand, lty=1, lwd=2, col="gray")
lines(x, under.random$results.pa, lty=1, lwd=2, col="black")
lines(x, under.random$results.obs, lty=6, lwd=2, col="black")

plot(c(0, 1), c(0, 1), type="n", xlab="Proportion of nodes removed", # Plot results
     ylab="Relative size of giant component",  main="(b) Overreporting")
lines(x, over.random$results.rand, lty=1, lwd=2, col="gray")
lines(x, over.random$results.pa, lty=1, lwd=2, col="black")
lines(x, over.random$results.obs, lty=6, lwd=2, col="black")

plot(c(0, 1), c(0, 1), type="n", xlab="Proportion of nodes removed", # Plot results
     ylab="Relative size of giant component",  main="(c) Desire-to-link")
lines(x, desire.random$results.rand, lty=1, lwd=2, col="gray")
lines(x, desire.random$results.pa, lty=1, lwd=2, col="black")
lines(x, desire.random$results.obs, lty=6, lwd=2, col="black")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 2, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Random", "PA", "Nyakatoke"), xpd = TRUE, horiz = TRUE, 
       inset = c(0, 0.02), lty=c(1, 1, 6), lwd=c(2, 2, 2), 
       col = c("gray", "black", "black"), cex = 1, text.width = c(0.2, 0.1, 0.2),
       seg.len=5)

# Robustness to targeted attack
m <- matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4), nrow=10, ncol=1)
layout(m)
par(mar = c(4, 4, 3, 1))

plot(c(0, 1), c(0, 1), type="n", xlab="Proportion of nodes removed", # Plot results
     ylab="Relative size of giant component", main="(a) Underreporting")
lines(x, under.targeted$results.rand, lty=1, lwd=2, col="gray")
lines(x, under.targeted$results.pa, lty=1, lwd=2, col="black")
lines(x, under.targeted$results.obs, lty=6, lwd=2, col="black")

plot(c(0, 1), c(0, 1), type="n", xlab="Proportion of nodes removed", # Plot results
     ylab="Relative size of giant component",  main="(b) Overreporting")
lines(x, over.targeted$results.rand, lty=1, lwd=2, col="gray")
lines(x, over.targeted$results.pa, lty=1, lwd=2, col="black")
lines(x, over.targeted$results.obs, lty=6, lwd=2, col="black")

plot(c(0, 1), c(0, 1), type="n", xlab="Proportion of nodes removed", # Plot results
     ylab="Relative size of giant component",  main="(c) Desire-to-link")
lines(x, desire.targeted$results.rand, lty=1, lwd=2, col="gray")
lines(x, desire.targeted$results.pa, lty=1, lwd=2, col="black")
lines(x, desire.targeted$results.obs, lty=6, lwd=2, col="black")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 2, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Random", "PA", "Nyakatoke"), xpd = TRUE, horiz = TRUE, 
       inset = c(0, 0.02), lty=c(1, 1, 6), lwd=c(2, 2, 2), 
       col = c("gray", "black", "black"), cex = 1, text.width = c(0.2, 0.1, 0.2),
       seg.len=5)

