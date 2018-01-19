#### Functions for assortativity and clustering simulations for "The structure 
#### of risk-sharing networks" by H. Henderson and A. Alam.
#### Note: The algorithm used in these simulations is discussed in detail in
#### Holme and Zhao's (2007) "Exploring the assortavity-clustering space of
#### a network's degree sequence."

# Find maximum assortativity coefficient
extreme.r <- function(network, max, stop){
  stopifnot(is.simple(network), is.logical(max), is.numeric(stop)) # Check inputs
  ntw1 <- network
  r1 <- igraph::assortativity_degree(ntw1) # Calculate initial assortativity coeff.
  ctr <- 1  # Initialize counter
  max_min <- ifelse(max, 1, -1)  # Used in if statement below to "toggle" condition
  while (ctr <= stop){
    ntw2 <- igraph::rewire(ntw1, with=keeping_degseq(loops=FALSE, niter=1)) # One trial
    r2 <- igraph::assortativity_degree(ntw2) # Calculate new assortativity coeff.
    if (r2*max_min > r1*max_min){   # If new value is an improvement ...
      r1 <- r2   # ... update values ...
      ntw1 <- ntw2
      ctr <- 1  # ... and reset counter
    } else {      # Otherwise keep current values ...
      ctr <- ctr + 1   # and increment counter
    }
  }
  return(r1)   # Return max. or min. coefficient
}

# Find min. or max. clustering coefficient given an assortativity interval
extreme.c <- function(network, max, lower, upper, stop){
  ntw1 <- network
  r1 <- igraph::assortativity_degree(ntw1) # Calculate initial assortativity coeff.
  c1 <- igraph::transitivity(ntw1, type="global", isolates="zero") # Calculate initial clustering coeff.
  ctr <- 1
  max_min <- ifelse(max, 1, -1)  # Used in if statement below to "toggle" condition
  while(ctr <= stop){
    ntw2 <- igraph::rewire(ntw1, with=keeping_degseq(loops=FALSE, niter=1)) # One trial
    r2 <- igraph::assortativity_degree(ntw2) # Calculate new assortativity coeff.
    c2 <- igraph::transitivity(ntw2, type="global", isolates="zero") # Calculate new clustering coeff.
    cond1 <- (r2 < lower && r2 > r1)  # Conditions for if statements
    cond2 <- (r2 >= upper && r2 < r1)
    cond3 <- (lower <= r2 && r2 < upper)
    cond4 <- c2*max_min > c1*max_min
    if (cond1 || cond2){  # If outside interval and moving in right direction ...
      r1 <- r2    # ... then update ...
      c1 <- c2
      ntw1 <- ntw2
    }
    if (cond3 && cond4){  # If inside interval and moving in right direction ...
      r1 <- r2     # then update ...
      c1 <- c2
      ntw1 <- ntw2
      ctr <- 1     # ... and reset counter
    }
    if (cond3 && !cond4){ # If inside interval and clustering not improving ...
      ctr <- ctr + 1     # ... then increment counter
    }
  }
  return(c1)  # Return max. or min. clustering coefficient
}

# Execute extreme.c repeatedly
repeat.c <- function(network, max, lower, upper, stop, reps){
  stopifnot(is.simple(network), is.logical(max), 
            is.numeric(c(lower, upper, stop, reps))) # Check inputs
  results <- sapply(1:reps, function(x) extreme.c(network, max, lower, upper, stop))
  overall.c <- ifelse(max, max(results), min(results)) 
  return(overall.c)
}

# Randomly "walk" in a given pixel for more uniform sampling
walk <- function(network, r.lower, r.upper, c.lower, c.upper, stop){
  ntw1 <- network
  ctr <- 1
  while (ctr <= stop){
    ntw2 <- igraph::rewire(ntw1, with=keeping_degseq(loops=FALSE, niter=1)) # One trial
    r2 <- igraph::assortativity_degree(ntw2) # Calculate new assortativity coeff.
    c2 <- igraph::transitivity(ntw2, type="global", isolates="zero") # Calculate new clustering coeff.
    cond1 <- (r2 >= r.lower && r2 < r.upper)
    cond2 <- (c2 >= c.lower && c2 < c.upper)
    if (cond1 && cond2){ # If inside interval ...
      ntw1 <- ntw2   # ... then update
      ctr <- ctr + 1
    }
  }
  return(ntw2)
}

# Generate data frame w/ valid pixels only
valid.pixels <- function(min.r, max.r, vmin.c, vmax.c, L){
  stopifnot(is.numeric(c(min.r, max.r, vmin.c, vmax.c, L))) # Check inputs
  stopifnot(length(vmin.c)==L,length(vmax.c)==L) # Check dimensions of inputs
  min.c <- min(vmin.c) # Get overall min. and max. clustering coeff.
  max.c <- max(vmax.c)
  cuts.r <- seq(min.r, max.r, length.out=L + 1) # List cutoff points for valid region
  cuts.c <- seq(min.c, max.c, length.out=L + 1)
  # Create vectors w/ lower and upper bounds for all pixels in valid region
  lpixel.r <- cuts.r[1 : L] 
  upixel.r <- cuts.r[2 : (L + 1)]
  lpixel.c <- cuts.c[1 : L]
  upixel.c <- cuts.c[2 : (L + 1)]
  # Create data frames where each row represents an interval in the valid region
  # Note: vmin.c and vmax.c are used below to locate the valid pixels
  valid.r <- data.table(lpixel.r, upixel.r, vmin.c, vmax.c) 
  valid.c <- data.table(lpixel.c, upixel.c)
  # Create database where each row represents a pixel with r and c bounds
  pixels.r <- valid.r[rep(1:L, each=L),]
  pixels.c <- valid.c[rep(1:L, times=L),]
  all.pixels <- cbind(pixels.r, pixels.c)
  # Locate valid pixels and drop vmin.c and vmax.c columns                  
  valid.pixels <- subset(all.pixels, upixel.c > vmin.c & lpixel.c < vmax.c)
  valid.pixels <- subset(valid.pixels, select=-c(vmin.c, vmax.c))
  # Also return all pixels without vmin.c and vmax.c
  all.pixels <- subset(all.pixels, select=-c(vmin.c, vmax.c))
  return(list(all.pixels, valid.pixels))
}

# Sample all pixels and calculate statistics for each
sampling <- function(samp, name, network, pixels, min.r, max.r, min.c, max.c, stop){
  #stopifnot(is.character(name), is.simple(network), dim(pixels)[2]==4,
  #          is.numeric(c(samp, min.r, max.r, min.c, max.c, stop))) # Check inputs
  pixels <- pixels[sample(1:nrow(pixels)),] # Randomize pixel order
  pixels[, "dist"] <- numeric()  # Add columns to store results
  pixels[, "component"] <- numeric()
  pixels[, "random.50"] <- numeric()
  pixels[, "targeted.50"] <- numeric()
  pixels[, "targeted.75"] <- numeric()
  #pixels <- as.matrix(pixels)  
  ntw1 <- network
  r1 <- igraph::assortativity_degree(ntw1) # Calculate initial assortativity coeff.
  c1 <- igraph::transitivity(ntw1, type="global", isolates="zero") # Calculate initial clustering coeff.
  for (octr in 1:nrow(pixels)){  # Run for each pixel
    lb.r <- pixels[octr, 1]  # Unpack lower and upper bounds for pixel
    ub.r <- pixels[octr, 2]
    lb.c <- pixels[octr, 3]
    ub.c <- pixels[octr, 4]
    center.r <- (lb.r + ub.r)/2  # Calculate pixel center
    center.c <- (lb.c + ub.c)/2
    rd1 <- ((r1 - center.r)/(max.r - min.r))^2 # Calculate components of distance formula
    cd1 <- ((c1 - center.c)/(max.c - min.c))^2
    d1 <- sqrt(rd1 + cd1)   # Calculate initital distance
    cond1 <- (r1 < lb.r) || (r1 >= ub.r)
    cond2 <- (c1 < lb.c) || (c1 >= ub.c)
    ictr <- 1   # Initialize inner counter
    while ((cond1 || cond2) && ictr <= 1e+05){  # While outside pixel ...
      ntw2 <- igraph::rewire(ntw1, with=keeping_degseq(loops=FALSE, niter=1)) # One trial
      r2 <- igraph::assortativity_degree(ntw2) # Calculate new assortativity coeff.
      c2 <- igraph::transitivity(ntw2, type="global", isolates="zero") # Calculate new clustering coeff.
      rd2 <- ((r2 - center.r)/(max.r - min.r))^2
      cd2 <- ((c2 - center.c)/(max.c - min.c))^2
      d2 <- sqrt(rd2 + cd2)   # Calculate new distance
      if (d2 < d1) {  # If moving in the right direction ...
        ntw1 <- ntw2  # ... update values
        r1 <- r2
        c1 <- c2
        d1 <- d2
        cond1 <- (r1 < lb.r) || (r1 >= ub.r) # Check if still outside pixel
        cond2 <- (c1 < lb.c) || (c1 >= ub.c)
      } 
      ictr <- ictr + 1
    }  
    if (ictr <= 1e+05){   # If exit was successful ...
      ntw3 <- walk(ntw1, lb.r, ub.r, lb.c, ub.c, stop) # Execute walk
      dist <- igraph::mean_distance(ntw3)  # Calculate stats
      component <- max(igraph::components(ntw3)$csize)   
      random.50 <- f.random(ntw3, 0.50, 1)
      targeted.50 <- f.targeted(ntw3, 0.50)
      targeted.75 <- f.targeted(ntw3, 0.75)
      pixels[octr, "dist"] <- dist # Store stats
      pixels[octr, "component"] <- component  
      pixels[octr, "random.50"] <- random.50
      pixels[octr, "targeted.50"] <- targeted.50
      pixels[octr, "targeted.75"] <- targeted.75
    } 
  }
  write.csv(pixels, file=paste0("./out/", name, samp, ".csv"), row.names=FALSE)
  gc()
  return(NULL)
}

# Calculate f-robustness of a network where removal is random
f.random <- function(network, f, reps){
  #stopifnot(is.simple(network), f>0, f<1, is.numeric(reps)) # Check inputs
  ntw1 <- network
  n <- length(igraph::V(ntw1))
  results <- sapply(1:reps, function(i){
    sapply(1:(n - 2), function(j){
      to.remove <- sample(1:n, j)
      ntw2 <- igraph::delete_vertices(ntw1, to.remove) 
      max(igraph::components(ntw2)$csize) # Store comp. size  
    })
  })
  results <- rowMeans(results) # Average across repetitions
  results <- results/max(igraph::components(ntw1)$csize) # Relative comp. size
  result <- min(which(results < f), 119)/n  # Fraction to delete
  return(result)
}  

# Calculate f-robustness of a network where removal is targeted
f.targeted <- function(network, f){
  #stopifnot(is.simple(network), f>0, f<1) # Check inputs
  ntw1 <- network
  n <- length(igraph::V(ntw1))
  mat <- matrix(ncol=2, nrow=n, 0)  # Create matrix w/ node degrees
  mat[, 1] <- 1:n
  deg <- igraph::degree(ntw1)
  mat[, 2] <- deg
  mat <- mat[order(mat[, 2], decreasing=TRUE), ]  # Order matrix by degree
  results <- sapply(1:(n - 2), function(i){
    ntw2 <- igraph::delete_vertices(ntw1, mat[1:i, 1])
    max(igraph::components(ntw2)$csize)
  })
  results <- results/max(igraph::components(ntw1)$csize) # Relative comp. size
  result <- min(which(results < f), 119)/n  # Fraction to delete
  return(result)
}







