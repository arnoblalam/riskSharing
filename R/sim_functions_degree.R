#### Functions for degree distrubtion simulations for "The structure of 
#### risk-sharing networks" by H. Henderson and A. Alam.
#### Note: The network generating process used in these simulations is based on 
#### Jackson and Rogers' (2007) "Meeting stragers and friends of friends: How 
#### random are social networks?"

# Generate growing random networks w/ given number of nodes and links
jr <- function(n, links, pr, pn, directed){
  #stopifnot(is.numeric(c(n, links, pr, pn)), is.logical(directed))
  m <- ceiling(links/n)  # mr = mn; ceiling ensures first entrant can be accommodated
  init <- 2*m + 1  # Min. size of network at initialization
  ntw <- igraph::make_full_graph(init, directed=directed)  # Fully connected at initialization
  resid.nodes <- n - init
  if (identical(pr, 0) & directed==TRUE){
    resid.links <- links - igraph::ecount(ntw) - resid.nodes  # Accomodate d0 links
  } else {
    resid.links <- links - igraph::ecount(ntw)
  }
  resid.seq <- deg.seq(resid.nodes, resid.links)
  for (i in 1:resid.nodes){
    ntw <- igraph::add_vertices(ntw, 1)
    if (identical(pr, 0) & directed==TRUE){ # Add d0 link
      d0 <- sample(1:(igraph::vcount(ntw) - 1), 1)
      ntw <- igraph::add_edges(ntw, c(d0, init + i))
    } 
    r.candidates <- sample(1:(igraph::vcount(ntw) - 1), resid.seq[i]) # Raises error if mr > vcount(ntw)
    union <- igraph::ego(ntw, order=1, r.candidates, mode="out") # Mode ignored for undirected networks
    union <- unique(unlist(union))
    union <- union[!union %in% r.candidates] # Delete parents from union
    n.candidates <- sample(union, min(resid.seq[i], length(union)))  
    if (length(union) < resid.seq[i]){
      resid.seq[i + 1] <- resid.seq[i + 1] + (resid.seq[i] - length(union)) # Add shortfall to next iteration
    }
    r.selected <- r.candidates[sample(c(TRUE, FALSE), length(r.candidates), replace=TRUE, prob=c(pr, 1 - pr))]
    n.selected <- n.candidates[sample(c(TRUE, FALSE), length(n.candidates), replace=TRUE, prob=c(pn, 1 - pn))]
    to.link <- unique(c(r.selected, n.selected))
    rep <- rep(init + i, length(to.link))
    to.link <- c(rbind(rep, to.link))
    ntw <- igraph::add_edges(ntw, to.link)
  }
  return(ntw)
}

# Generate uniform degree sequence given the number of nodes and links
deg.seq <- function(n, links){
  #stopifnot(is.numeric(c(n, links)))  # Check inputs
  vec <- rep(round(links/n), n)  # Create uniform vector
  diff <- links - sum(vec)    
  samp <- sample(n, abs(diff))  # Add/subtract from random vector elements
  for (i in samp){
    vec[i] <- vec[i] + sign(diff)
  }
  if (any(vec < 0)) stop('Nodes cannot enter with negative number of links')
  if (sum(vec) != links) stop('Out sequence is incorrect')
  return(vec)
}

# Robustness to random attack
random <- function(network, reps){
  #stopifnot(is.igraph(network), is.numeric(reps))  # Check inputs
  ntw1 <- network
  n <- igraph::gorder(ntw1) 
  comp1 <- max(igraph::components(ntw1)$csize)
  results <- sapply(1:reps, function(i){
    sapply(1:(n - 2), function(j){
      to.remove <- sample(1:n, j)
      ntw2 <- igraph::delete_vertices(ntw1, to.remove) 
      max(igraph::components(ntw2)$csize)/comp1 # Calculate relative component size  
    })
  })
  results <- rowMeans(results) # Average across repetitions
  return(results)  
}
  
# Robustness to targeted attack
targeted <- function(network, directed){
  #stopifnot(is.igraph(network), is.logical(directed))  # Check inputs
  ntw1 <- network
  n <- igraph::gorder(ntw1)
  comp1 <- max(igraph::components(ntw1)$csize)
  mat <- matrix(ncol=2, nrow=n, 0)  # Create matrix w/ node degrees
  mat[, 1] <- 1:n
  if (directed==TRUE){
    mat[, 2] <- igraph::degree(ntw1, mode="in")
  } else {
    mat[, 2] <- igraph::degree(ntw1)  
  }
  mat <- mat[order(mat[, 2], decreasing=TRUE), ]  # Order matrix by degree
  results <- sapply(1:(n - 2), function(i){
    ntw2 <- igraph::delete_vertices(ntw1, mat[1:i, 1])
    max(igraph::components(ntw2)$csize)/comp1
  })
  return(results)
}

# Function for running simulations
degree.sims <- function(network, attack, reps){
  stopifnot(is.igraph(network), is.numeric(reps), 
            attack == "targeted" | attack =="random") # Check inputs
  ntw.obs <- network
  n <- igraph::vcount(ntw.obs) # Network characteristics
  links <- igraph::ecount(ntw.obs)
  directed <- igraph::is_directed(ntw.obs)
  results.obs <- matrix(ncol=(n - 2), nrow=reps, 0)  # Matrices to store results
  results.rand <- matrix(ncol=(n - 2), nrow=reps, 0) 
  results.pa <- matrix(ncol=(n - 2), nrow=reps, 0) 
  for (i in 1:reps){
    # Generate random network
    ntw.rand <- jr(n, links, pr=1, pn=0, directed)
    # Generate preferential attachment network
    ntw.pa <- jr(n, links, pr=0, pn=1, directed)
    if (attack=="random"){    # Generate results
      results.obs[i,] <- random(ntw.obs, 1)
      results.rand[i,] <- random(ntw.rand, 1)
      results.pa[i,] <- random(ntw.pa, 1)
    } 
    if (attack=="targeted"){
      results.obs[i,] <- targeted(ntw.obs, directed)
      results.rand[i,] <- targeted(ntw.rand, directed)
      results.pa[i,] <- targeted(ntw.pa, directed)
    }
  }
  results.obs <- colMeans(results.obs)  # Average across repetitions
  results.rand <- colMeans(results.rand)
  results.pa <- colMeans(results.pa)
  results <- data.frame(results.obs, results.rand, results.pa)
  return(results)
}
