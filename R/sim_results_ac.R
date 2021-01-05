#### Project title: The structure of risk-sharing networks
#### Authors: H. Henderson and A. Alam
#### Last updated: 1/1/2021
#### Purpose: Creates graphs for Section 4.2 of the paper

# Set up and read in data
# setwd("~/Documents/Articles/Risk-Sharing-Networks/Code")
library(igraph)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
all.over <- read.csv("./out_final/all_overreporting.csv")
all.under <- read.csv("./out_final/all_underreporting.csv")
nyakatoke <- read.csv("Nyakatoke.csv")
set.seed(12345)

# Load overreporting results and combine into one data frame
over.files <- list.files(path="./out_final", pattern="overreporting", full.names=TRUE)
over.files <- over.files[over.files != "./out_final/all_overreporting.csv"]
sims.over <- read.csv("./out_final/overreporting1.csv")[, 1:4]
ctr <- 1
for (i in over.files){
  draw <- read.csv(i)
  colnames(draw)[5] <- paste0("dist", ctr) # Rename columns for merge
  colnames(draw)[6] <- paste0("component", ctr)
  colnames(draw)[7] <- paste0("random.50", ctr)
  colnames(draw)[8] <- paste0("targeted.50", ctr)
  colnames(draw)[9] <- paste0("targeted.75", ctr)
  sims.over <- merge(sims.over, draw, by=c("lpixel.r","upixel.r", "lpixel.c", "upixel.c"))
  ctr <- ctr + 1
}
colnames=c("lpixel.r", "upixel.r", "lpixel.c", "upixel.c", "dist", "component",
           "random.50", "targeted.50", "targeted.75")
sims.over <- sapply(colnames, function(x) rowMeans(sims.over[,   # Calculate row means
                              grep(x, names(sims.over)), drop=FALSE], na.rm=TRUE))

# Repeat for underreporting network
under.files <- list.files(path="./out_final", pattern="underreporting", full.names=TRUE)
under.files <- under.files[under.files != "./out_final/all_underreporting.csv"]
sims.under <- read.csv("./out_final/underreporting1.csv")[, 1:4]
ctr <- 1
for (i in under.files){
  draw <- read.csv(i)
  colnames(draw)[5] <- paste0("dist", ctr) # Rename columns for merge
  colnames(draw)[6] <- paste0("component", ctr)
  colnames(draw)[7] <- paste0("random.50", ctr)
  colnames(draw)[8] <- paste0("targeted.50", ctr)
  colnames(draw)[9] <- paste0("targeted.75", ctr)
  sims.under <- merge(sims.under, draw, by=c("lpixel.r","upixel.r", "lpixel.c", "upixel.c"))
  ctr <- ctr + 1
}
sims.under <- sapply(colnames, function(x) rowMeans(sims.under[,   # Calculate row means
                            grep(x, names(sims.under)), drop=FALSE], na.rm=TRUE))

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

# Fix missing node issue
missing.vertices <- c(7, 30, 32, 46, 44, 65, 84, 88, 91, 96, 107, 110,
                      116, 117, 118, 119)
g.overreporting <- (Reduce(f = function(x, y) {y + igraph::vertex(x)},
                           x = missing.vertices,
                           init = g.overreporting,
                           right = TRUE))

# Remove multiple edges from overreporting network
g.overreporting <- igraph::simplify(g.overreporting)

# Select network statistics
r.under <- igraph::assortativity_degree(g.underreporting)
r.over <- igraph::assortativity_degree(g.overreporting)
c.under <- igraph::transitivity(g.underreporting, type="global")
c.over <- igraph::transitivity(g.overreporting, type="global")

# Merge dataframes
df.over <- merge(all.over, sims.over, by=c("lpixel.r","upixel.r", "lpixel.c",
                                          "upixel.c"), all.x=TRUE)
df.over["component"] <- df.over["component"]/119  # Normalize component
df.under <- merge(all.under, sims.under, by=c("lpixel.r","upixel.r", "lpixel.c",
                                          "upixel.c"), all.x=TRUE)
df.under["component"] <- df.under["component"]/119 # Normalize component

# Manipulate dataframes
L <- sqrt(length(df.over[, 1]))
df.over <- df.over[with(df.over, order(lpixel.r, lpixel.c)),]  # Sort dataframes
df.under <- df.under[with(df.under, order(lpixel.r, lpixel.c)),]
df.over$r <- rep(1:L, each=L)  # Assign each r and c interval a number 1:L ...
df.under$r <- rep(1:L, each=L)  # ... this is helpful for labeling the heatmaps.
df.over$c <- rep(1:L, times=L)
df.under$c <- rep(1:L, times=L)

# Retrieve min. and max. r and c values
min.r.over <- min(df.over$lpixel.r)
max.r.over <- max(df.over$upixel.r)
min.c.over <- min(df.over$lpixel.c)
max.c.over <- max(df.over$upixel.c)
min.r.under <- min(df.under$lpixel.r)
max.r.under <- max(df.under$upixel.r)
min.c.under <- min(df.under$lpixel.c)
max.c.under <- max(df.under$upixel.c)

# Get stats from random sample of r-c space for underreporting network
draws <- 1e+04
chains <- 10
r.results <- matrix(0, nrow=chains, ncol=draws)
c.results <- matrix(0, nrow=chains, ncol=draws)
for (i in 1:chains){
    print(paste("Entering chain", i))
    # Start each chain w/ a different randomly rewired network
    ntw1 <- igraph::rewire(g.underreporting, with=keeping_degseq(loops=FALSE, niter=1000))
    for (j in 1:draws){
        ntw2 <- igraph::rewire(ntw1, with=keeping_degseq(loops=FALSE, niter=1)) # One trial
        r.results[i, j] <- igraph::assortativity_degree(ntw2)
        c.results[i, j] <- igraph::transitivity(ntw2, type="global")
        ntw1 <- ntw2
    }
}
r.under.avg <- mean(r.results)
r.under.std <- sd(r.results)
c.under.avg <- mean(c.results)
c.under.std <- sd(c.results)

# Get stats from random sample of r-c space for overreporting network
r.results <- matrix(0, nrow=chains, ncol=draws)
c.results <- matrix(0, nrow=chains, ncol=draws)
for (i in 1:chains){
  print(paste("Entering chain", i))
  # Start each chain w/ a different randomly rewired network
  ntw1 <- igraph::rewire(g.overreporting, with=keeping_degseq(loops=FALSE, niter=1000))
  for (j in 1:draws){
    ntw2 <- igraph::rewire(ntw1, with=keeping_degseq(loops=FALSE, niter=1)) # One trial
    r.results[i, j] <- igraph::assortativity_degree(ntw2)
    c.results[i, j] <- igraph::transitivity(ntw2, type="global")
    ntw1 <- ntw2
  }
}
r.over.avg <- mean(r.results)
print(r.over.avg)
r.over.std <- sd(r.results)
print(r.over.std)
c.over.avg <- mean(c.results)
print(c.over.avg)
c.over.std <- sd(c.results)
print(c.over.std)
(r.over - r.over.avg)/r.over.std
(c.over - c.over.avg)/c.over.std

# Create axis labels
at.x <- seq(0.5, L + 0.5, length.out=5) # Must be offset to label at endpoints
at.y <- seq(0.5, L + 0.5, length.out=5)
label.x.over <- sprintf("%.2f", round(seq(min.r.over, max.r.over, length.out=5), 2))
label.y.over <- sprintf("%.2f", round(seq(min.c.over, max.c.over, length.out=5), 2))
label.x.under <- sprintf("%.2f", round(seq(min.r.under, max.r.under, length.out=5), 2))
label.y.under <- sprintf("%.2f", round(seq(min.c.under, max.c.under, length.out=5), 2))

# Coordinates in r-c space for both the actual and the expected network
r.aloc.under <- 0.5 + L*(r.under - min.r.under)/(max.r.under - min.r.under)
c.aloc.under <- 0.5 + L*(c.under - min.c.under)/(max.c.under - min.c.under)
r.aloc.over <- 0.5 + L*(r.over - min.r.over)/(max.r.over - min.r.over)
c.aloc.over <- 0.5 + L*(c.over - min.c.over)/(max.c.over - min.c.over)

r.eloc.under <- 0.5 + L*(r.under.avg - min.r.under)/(max.r.under - min.r.under)
c.eloc.under <- 0.5 + L*(c.under.avg - min.c.under)/(max.c.under - min.c.under)
r.eloc.over <- 0.5 + L*(r.over.avg - min.r.over)/(max.r.over - min.r.over)
c.eloc.over <- 0.5 + L*(c.over.avg - min.c.over)/(max.c.over - min.c.over)

# Create 4x4 plot for overreporting results
plot1 <- levelplot(dist ~ r*c,
    df.over,
    col.regions=gray(80:0/100),
    main="(a) Average Distance",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.over),
    y=list(at=at.y, label=label.y.over))) +
xyplot(c.aloc.over ~ r.aloc.over, pch=4, col="white") +
xyplot(c.eloc.over ~ r.eloc.over, pch=3, col="white")
plot2 <- levelplot(component ~ r*c,
    df.over,
    col.regions=gray(95:0/100),
    main="(b) Size of Giant Component",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.over),
    y=list(at=at.y, label=label.y.over))) +
xyplot(c.aloc.over ~ r.aloc.over, pch=4, col="white") +
xyplot(c.eloc.over ~ r.eloc.over, pch=3, col="white")
plot3 <- levelplot(random.50 ~ r*c,
    df.over,
    col.regions=gray(80:0/100),
    main="(c) Robustness to Random Failure",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.over),
    y=list(at=at.y, label=label.y.over))) +
xyplot(c.aloc.over ~ r.aloc.over, pch=4, col="white") +
xyplot(c.eloc.over ~ r.eloc.over, pch=3, col="white")
plot4 <- levelplot(targeted.50 ~ r*c,      # Overreporting network
    df.over,
    col.regions=gray(80:0/100),
    main="(d) Robustness to Targeted Attack",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.over),
    y=list(at=at.y, label=label.y.over))) +
xyplot(c.aloc.over ~ r.aloc.over, pch=4, col="white") +
xyplot(c.eloc.over ~ r.eloc.over, pch=3, col="white")
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)

# Create 4x4 plot for underreporting results
plot1 <- levelplot(dist ~ r*c,
    df.under,
    col.regions=gray(80:0/100),
    main="(a) Average Distance",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.under),
    y=list(at=at.y, label=label.y.under))) +
xyplot(c.aloc.under ~ r.aloc.under, pch=4, col="white") +
xyplot(c.eloc.under ~ r.eloc.under, pch=3, col="white")
plot2 <- levelplot(component ~ r*c,
    df.under,
    col.regions=gray(95:0/100),
    main="(b) Size of Giant Component",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.under),
    y=list(at=at.y, label=label.y.under))) +
xyplot(c.aloc.under ~ r.aloc.under, pch=4, col="white") +
xyplot(c.eloc.under ~ r.eloc.under, pch=3, col="white")
plot3 <- levelplot(random.50 ~ r*c,
    df.under,
    col.regions=gray(80:0/100),
    main="(c) Robustness to Random Failure",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.under),
    y=list(at=at.y, label=label.y.under))) +
xyplot(c.aloc.under ~ r.aloc.under, pch=4, col="white") +
xyplot(c.eloc.under ~ r.eloc.under, pch=3, col="white")
plot4 <- levelplot(targeted.50 ~ r*c,
    df.under,
    col.regions=gray(80:0/100),
    main="(d) Robustness to Targeted Attack",
    xlab="Assortativity",
    ylab="Clustering",
    scale=list(
    x=list(at=at.x, label=label.x.under),
    y=list(at=at.y, label=label.y.under))) +
xyplot(c.aloc.under ~ r.aloc.under, pch=4, col="white") +
xyplot(c.eloc.under ~ r.eloc.under, pch=3, col="white")
grid.arrange(plot1, plot2, plot3, plot4, ncol=2)


