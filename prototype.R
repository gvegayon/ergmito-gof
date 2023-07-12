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
  seq2 <- sort(unique(S.$statmat[,2]))
  
  # Computing the conditional probability
  ans <- vector("list", length(seq2))
  for (i in seq_along(ans))
    ans[[i]] <- d_ergm(a. = seq1, b. = seq2[i], S. = S., theta. = theta)
  
  # Calculating the quantiles
  quantiles <- sapply(ans, function(a) {
    
    a <- cumsum(a)
    
    # Identifying which are zero
    are0 <- which(a < 0.000000001)
    
    if (length(are0)) {
      a <- a[-are0]
      e <- seq1[-are0]
    } else
      e <- seq1
    
    if (!length(e))
      return(c(NA, NA, NA))
    
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
  
  # Which quantile of the conditional is closer to the 50%
  ergm_probs <- ergmito::exact_loglik(
    S.$statmat,
    params        = theta,
    stats_weights = S.$weights,
    stats_statmat = S.$statmat,
    as_prob       = TRUE
    ) * S.$weights
  
  ord <- order(S.$statmat[,2], decreasing = FALSE)
  ergms_probs <- ergm_probs[ord]
  
  b_50pcent <- which.min(abs(cumsum(ergm_probs) - .5))
  b_50pcent <- S.$statmat[ord,][b_50pcent,2]
  b_max     <- S.$statmat[ord,][which.max(ergms_probs),2]
  
  structure(
    list(
      a = seq1,
      b = seq2,
      pr = ans,
      quantiles = cbind(
        `5%`= quantiles[,1],
        `50%` = quantiles[,2],
        `95%` = quantiles[,3]
        ),
      par_names = colnames(S.$statmat),
      probs2    = ergm_probs,
      b_50pcent = b_50pcent,
      b_max     = b_max,
      support   = S.
      ),
    class = "ergm_conditional"
  )
}

namer <- function(x) {
  switch(
    x, 
    edges = "Edge count",
    ttriple = "Transitive Triads",
    mutual = "Mutual Ties",
    nodeicov.gender = "Gender-Receiver Effect",
    nodematch.gender = "Gender-Homophily"
  )
}

#' Plot method for conditional probabilities on ERGMs
#' @noRd
plot.ergm_conditional <- function(x, xlab = namer(x$par_names[2]), ylab = namer(x$par_names[1]), ...) {
  
  # scaler <- function(a) {
  #   (a - min(x$b))/diff(range(x$b))
  # }
  scaler <- function(a) a
  
  y_range <- range(x$a)
  x_range <- scaler(range(x$b))
  
  # # Border
  # dat <- data.frame(
  #   x = c(x$b, rev(x$b)),
  #   y = c(x$quantile[,1], rev(x$quantile[,3]))
  # )
  # 
  # ggplot2::ggplot(data = dat, aes(x=x,y=y)) +
  #   ggplot2::geom_polygon()
  # 
  graphics::plot(
    NA,
    xlim = x_range,
    ylim = y_range,
    xlab = xlab,
    ylab = ylab
    )
  
  op <- graphics::par(xpd = FALSE)
  graphics::grid()
  graphics::par(op)
  
  # Drawing a polygon
  graphics::polygon(
    x = scaler(c(x$b, rev(x$b))),
    y = c(x$quantile[,1], rev(x$quantile[,3])),
    col = "lightgray",
    border = "darkgray"
      
  )
  
  lines(
    y    = x$quantiles[,2],
    x    = scaler(x$b),
    type = "l", col = "red", lwd = 2, lty = 2
  )
  
  # abline(v = x$b_50pcent, lty = 2, lwd = 2)
  
  return()
  
}

# Plot the first with a title
plot_first <- function(x, main="", ...) {
  for (i in seq_along(x)) {
    plot(x[[i]], ...)
    if (i == 1L) {
      op <- par(xpd=NA)
      # title(main = main, line = 1.25, font.main = 1, cex.main = 1.5)
      par(op)
    }
  }

}

# Testing the function
# Conditional probabilities
theta <- c(.5, -.5)

# Stats counts
x <- network(matrix(0, 5, 5))
set.vertex.attribute(x, "gender", c(0,1,0,1,0))


# General plotting parameters
parpar <- list(
  mfcol = c(2, 4), mar = c(4,2.5,.5,1.5), oma = c(0,2,2.25,0)
)
width. <- 7
height. <- 3.5


# Mutual --------------------------------------
models_mutual <- list(
  x ~ mutual + edges,
  x ~ mutual + ttriad,
  x ~ mutual + nodematch("gender"),
  x ~ mutual + nodeicov("gender")
)

ans_same <- lapply(models_mutual, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_mutual, conditional_dist, theta = c(2, 0))

graphics.off()
svg("conditional-prob-mutuals.svg", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = "(a)")
plot_first(ans_diff, main = "(b)")
par(op)
title(ylab = "Number of Mutual ties")
dev.off()

# Triads --------------------------------------
models_ttriad <- list(
  x ~ ttriad + edges,
  x ~ ttriad + mutual,
  x ~ ttriad + nodematch("gender"),
  x ~ ttriad + nodeicov("gender")
)

ans_same <- lapply(models_ttriad, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_ttriad, conditional_dist, theta = c(1, 0))

graphics.off()
svg("conditional-prob-ttriad.svg", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = "(a)")
plot_first(ans_diff, main = "(b)")
par(op)
title(ylab = "Number of Transitive Triads")
dev.off()

# Gender homophily --------------------------------------
models_homophily <- list(
  x ~ nodematch("gender") + edges,
  x ~ nodematch("gender") + mutual,
  x ~ nodematch("gender") + ttriad,
  x ~ nodematch("gender") + nodeicov("gender")
)

ans_same <- lapply(models_homophily, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_homophily, conditional_dist, theta = c(2, 0))

graphics.off()
svg("conditional-prob-homophily.svg", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = "(a)")
plot_first(ans_diff, main = "(b)")
par(op)
title(ylab = "Number of Gender-Homophilic Ties")
dev.off()

# Gender homophily --------------------------------------
models_icov <- list(
  x ~ nodeicov("gender") + edges,
  x ~ nodeicov("gender") + mutual,
  x ~ nodeicov("gender") + ttriad,
  x ~ nodeicov("gender") + nodematch("gender")
)

ans_same <- lapply(models_icov, conditional_dist, theta = c(0, 0))
ans_diff <- lapply(models_icov, conditional_dist, theta = c(2, 0))

graphics.off()
svg("conditional-prob-receiver-effect.svg", width = width., height = height.)
op <- do.call(par, parpar)
plot_first(ans_same, main = "(a)")
plot_first(ans_diff, main = "(b)")
par(op)
title(ylab = "Gender-Receiver Effect")
dev.off()

