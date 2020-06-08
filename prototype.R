library(ergmito)
library(ergm)

#' Generates the conditional probability of observing e given (t, theta.)
#' @noRd
d_ergm <- function(a., b., S., theta.) {
  
  if (length(a.) > 1)
    return(sapply(a., d_ergm, b. = b., theta. = theta., S. = S.))
  
  # Which matches
  idx <- which(colSums(t(S.$statmat) == c(a., b.)) == 2)
  if (!length(idx))
    return(0)
  
  idx0 <- which(S.$statmat[,2] == b.)
  
  # Computing the probability
  (S.$weights[idx] * exp(theta. %*% cbind(c(a., b.)))/
    rbind(S.$weights[idx0]) %*% exp(S.$statmat[idx0, , drop=FALSE] %*% theta.))[1]
  
}

#' Generates the c(.05,.5,.95) quantiles for the conditional distribution
#' of term 1 given term 2.
#' @param model A model to be passed to [ergm::ermg.allstats]
#' @param theta Model parameters
conditional_dist <- function(model, theta) {
  S. <- ergm::ergm.allstats(model)
  
  seq1 <- sort(unique(S.$statmat[,1]))
  seq2 <- sort(unique(S.$statmat[,1]))
  
  # Computing the conditional probability
  ans <- vector("list", length(seq2))
  for (i in seq_along(ans))
    ans[[i]] <- d_ergm(a. = seq1, b. = seq2[i], S. = S., theta. = theta)
  
  # Calculating the quantiles
  quantiles <- sapply(ans, function(a) {
    
    a <- cumsum(a)
    
    # Identifying which are zero
    are0 <- which(a < 0.000000001)
    
    # if (length(are0)) {
    #   a <- a[-are0]
    #   e <- seq1[-are0]
    # } else
    #   e <- seq1
    
    e[c(
      # 5% 
      which.min(abs(a - .05)),
      # 50%
      which.min(abs(a - .5)),
      # %95
      which.min(abs(a - .95))
    )]
    
  })
  quantiles <- t(quantiles)
  
  list(
    x = seq2,
    y = cbind(
      `5%`= quantiles[,1],
      `50%` = quantiles[,2],
      `95%` = quantiles[,3]
      )
  )
}

# Testing the function
# Conditional probabilities
theta <- c(.5, -.5)

# Stats counts
x <- network(matrix(0, 5, 5))
set.vertex.attribute(x, "gender", c(0,1,0,1,0))
conditional_dist(x ~ mutual + nodematch("gender"), theta = theta)

# Getting the ranges
r_edges <- sort(unique(S$statmat[,1]))
r_ttriad <- range(S$statmat[,2])

r_ttriad <- sort(unique(S$statmat[,2]))
ans <- vector("list", length(r_ttriad))
for (i in seq_along(ans)) {
  ans[[i]] <- d_ergm(r_edges, r_ttriad[i])
}

# Plotting:  -------------------------------------------------------------------

# 5%
x <- sapply(ans, function(a) {
  
  a <- cumsum(a)
  
  # Identifying which are zero
  are0 <- which(a < 0.000000001)
  
  if (length(are0)) {
    a <- a[-are0]
    e <- r_edges[-are0]
  } else
    e <- r_edges

  e[c(
    # 5% 
    which.min(abs(a - .05)),
    # 50%
    which.min(abs(a - .5)),
    # %95
    which.min(abs(a - .95))
  )]

})

plot(x = r_ttriad, y = x[2,], type = "b", ylim = range(r_edges), lty = 2,
     xlab = colnames(S$statmat)[2], ylab = colnames(S$statmat)[1])
lines(x = r_ttriad, y = x[1,], type = "l", col = "red", lwd = 2, lty = 2)
lines(x = r_ttriad, y = x[3,], type = "l", col = "blue", lwd = 2, lty = 2)
grid(nx = length(r_ttriad), length(r_edges))
