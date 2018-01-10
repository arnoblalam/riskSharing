## ----load-data-----------------------------------------------------------
# Load package
library(igraph)
library(riskSharing)

#Get data
data("nyakatoke")

# Create the underreporting graph

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

## ----percolation---------------------------------------------------------
S <- numeric(120) # vector to hold our results
sim_size <- 1000 # Number of simulations
i <- numeric(sim_size) # Vector to hold cluster sizes
for (n in 0:119) {
    # n is the number of nodes to remove
    for (m in 1:sim_size) {
        # m is the mth simulation
        to_remove <- sample(1:119, n)
        # randomly pick n of the 119 nodes to remove
        i[m] <- components(
            delete_vertices(g.underreporting, to_remove))$csize[1]/119
        # calculate the size of the giant component and save the result to
        # the variable i
    }
    S[n+1] <- mean(i)
    # the average cluster size (after removing n nodes) is the average of the 
    # i variable
}

## ----phi v S underrepprting, fig.width=6, fig.height=4-------------------
my_df <- data.frame(phi = 1 - (0:119)/119, S = S)
# rbind(my_df, c(phi = 0, S = 1), c(phi = 1, S = 0))
main = expression(paste("Occupaion Probability, ", phi, ", vs. size of giant cluster ", S))
plot(my_df[1:117,], 
     type = "l", 
     main = main,
     xlab = expression(phi))
# We will also add the theoretical results - first exponential
# library(gsl) # For Lambert's W function
lambda <- mean(degree(g.underreporting))
plot(function(phi) phi+gsl::lambert_W0(-lambda * phi * exp(-lambda*phi))/lambda, from = 0, to = 1, add = TRUE, 
     col = "red")

# Now the power law theoretical results.  We first obtain the MLE estimates
# alpha <- power.law.fit(degree(g.underreporting))$alpha
# k <- degree(g.underreporting)
# Ek <- mean(k)
# Ek2 <- mean(k^2)
# phi_c <- Ek/(Ek2 - Ek)
# plot(function(phi) (2*Ek/0.810891)*(phi  - phi_c)/phi_c, from = 0, to = 1, add = TRUE, col = "blue")
legend("topleft", legend = c("Simulation", "Poisson"), col = c("black", "red"), lty = c(1, 1))

## ----overreporting model create graph------------------------------------
# create the graph (see the basic statistics section for  documentation of the
# procedure below)
overreporting.df <- 
    nyakatoke[nyakatoke$willingness_link1 == 1 & nyakatoke$willingness_link2 == 1,]
g.overreporting <- igraph::graph_from_data_frame(overreporting.df, directed = FALSE)
missing.vertices <- c(7, 30, 32, 46, 44, 65, 84, 88, 91, 96, 107, 110, 
                      116, 117, 118, 119)
g.overreporting <- (Reduce(f = function(x, y) {y + igraph::vertex(x)},
       x = missing.vertices,
       init = g.overreporting,
       right = TRUE))
g.overreporting <- igraph::simplify(g.overreporting)

## ----percolation overreporting, fig.width=6, fig.height=4----------------
S <- numeric(120) # vector to hold our results
sim_size <- 1000 # Number of simulations
i <- numeric(sim_size) # Vector to hold cluster sizes
for (n in 0:119) {
    # n is the number of nodes to remove
    for (m in 1:sim_size) {
        # m is the mth simulation
        to_remove <- sample(1:119, n)
        # randomly pick n of the 119 nodes to remove
        i[m] <- components(
            delete_vertices(g.overreporting, to_remove))$csize[1]/119
        # calculate the size of the giant component and save the result to
        # the variable i
    }
    S[n + 1] <- mean(i)
    # the average cluster size (after removing n nodes) is the average of the 
    # i variable
}
my_df <- data.frame(phi = 1 - (0:119)/119, S = S)
# rbind(my_df, c(phi = 0, S = 1), c(phi = 1, S = 0))
main = expression(paste("Occupaion Probability, ", phi, ", vs. size of giant cluster ", S))
plot(my_df[1:117,], 
     type = "l", 
     main = main,
     sub = "Overreporting Model",
     xlab = expression(phi))
lambda <- mean(degree(g.overreporting))
plot(function(phi) phi + gsl::lambert_W0(-lambda * phi * exp(-lambda*phi))/lambda, from = 0, to = 1, add = TRUE, 
     col = "red")


legend("topleft", legend = c("Simulation", "Poisson"), col = c("black", "red"), lty = c(1, 1))

