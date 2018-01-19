#### Assortativity and clustering simulations for "The structure of risk-sharing
#### networks" by H. Henderson and A. Alam.
#### Note: The algorithm used in these simulations is discussed in detail in
#### Holme and Zhao's (2007) "Exploring the assortavity-clustering space of
#### a network's degree sequence."

# Set up and read in data
# Note: igraph must be installed
setwd("/home/aa2288a")
library("parallel")
library("igraph") 
library("data.table", lib.loc="/home/aa2288a/R-packages/")
source("sim_functions_ac.R")
nyakatoke <- read.csv("Nyakatoke.csv")

# Underreporting network
underreporting.df <- 
  nyakatoke[nyakatoke$willingness_link1 == 1 | nyakatoke$willingness_link2 == 1,]

# Read underreporting network into igraph and remove redundant links
g.underreporting <- graph_from_data_frame(underreporting.df, directed = FALSE)
g.underreporting <- simplify(g.underreporting)

# Overreporting network
overreporting.df <- 
  nyakatoke[nyakatoke$willingness_link1 == 1 & nyakatoke$willingness_link2 == 1,]

# Read overreporting network into igraph
g.overreporting <- graph_from_data_frame(overreporting.df, directed = FALSE)

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
g.overreporting <- simplify(g.overreporting)

# Specify network to use
network <- g.underreporting  # Underreporting or overreporting 
name <- ifelse(identical_graphs(network, g.overreporting), "overreporting", 
               "underreporting")

# Create directory for results 
# Note: If the directory already exists, then delete results associated with 
# the network being analyzed.
dir.create("out", showWarnings=FALSE)     
file.remove(list.files(path="./out", pattern=name, full.names=TRUE))

# Remove select objects
rm(g.overreporting, g.underreporting, missing.vertices, nyakatoke,
   overreporting.df, underreporting.df)

# Parameters
# Note: Holme and Zhao parameter values are in parentheses
vsame <- 1e+02  # No. of times r and c don't change before stopping (10^5)
vrep <- 2   # No. of times extreme methods are called (5)
vsamp <- 2  # No. of times each pixel is sampled (100)
vrnd <- 10  # No. of discarded networks in walk method (1000)
L <- 5     # No. of partitions on both assortativity and clustering space (50) 
no.cores <- 20  # No. of cores to use

# Set seed
RNGkind("L'Ecuyer-CMRG") 
set.seed(1)

# Get min. and max. assortativity coefficients
min.r <- mclapply(1:vrep, function(x) extreme.r(network=network, max=FALSE, 
                                           stop=vsame), mc.cores = no.cores)
min.r <- min(unlist(min.r))  # Get overall minimum
print(min.r)
max.r <- mclapply(1:vrep, function(x) extreme.r(network=network, max=TRUE, 
                                           stop=vsame), mc.cores = no.cores)
max.r <- max(unlist(max.r))  # Get overall maximum
print(max.r)

# Get min. and max. clustering coefficients for each assortativity interval
cuts.r <- seq(min.r, max.r, length.out=L + 1) # List cutoff points for valid region
lpixel.r <- cuts.r[1 : L]   
upixel.r <- cuts.r[2 : (L + 1)]
valid.r <- data.frame(lpixel.r, upixel.r) # Data frame with assortativity intervals
vmin.c <- mclapply(1:L, function(x) repeat.c(network=network, max=FALSE, 
          lower=valid.r[x, 1], upper=valid.r[x, 2], stop=vsame, reps=vrep), 
          mc.cores = no.cores)
vmin.c <- unlist(vmin.c)
print(vmin.c)
vmax.c <- mclapply(1:L, function(x) repeat.c(network=network, max=TRUE, 
          lower=valid.r[x, 1], upper=valid.r[x, 2], stop=vsame, reps=vrep),
          mc.cores = no.cores)
vmax.c <- unlist(vmax.c)
print(vmax.c)

# Get valid pixels
pixels <- valid.pixels(min.r=min.r, max.r=max.r, vmin.c=vmin.c, vmax.c=vmax.c, L=L)
write.csv(pixels[[1]], file=paste0("./out/", "all_", name, ".csv"), row.names=FALSE)
pixels <- pixels[[2]] # Valid pixels

# Remove select objects
rm(cuts.r, extreme.c, extreme.r, L, lpixel.r, repeat.c, upixel.r, valid.pixels,
   valid.r, vrep, vsame)

# Repeatedly sample valid pixels
# Note: This procedure raises some errors from failed repetitions. 
min.c <- min(vmin.c) # Get overall min. and max. clustering coefficients
max.c <- max(vmax.c) 
seq <- seq(from=1, to=vsamp, by=no.cores)
for (i in seq){   # Explicitly enumerate iterations
  start <- i
  end <- i + no.cores - 1
  results <- mclapply(start:end, function(samp) sampling(samp, name, network, 
                pixels, min.r, max.r, min.c, max.c, vrnd), mc.cores = no.cores)
}
