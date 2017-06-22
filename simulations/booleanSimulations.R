source('./booleanUtilities.R')
source('./makeMat.R')
library(iRF)
n.cores <- 24

x.params <- list()
x.params$dist <- 'cauchy'
y.params <- list()
y.params$classification <- TRUE

rules <- 'xor' #one of: 'and', 'or', 'xor', 'mixxor'
thresholds <- 2
noises <- c(0.1, 0.15, 0.2)
p <- 100
n.obs <- 5000
active1 <- 1:8
active2 <- NULL
mixture.prop <- NULL

for (r in rules) {
  y.params$rule <- r
  for (thresh in thresholds) {
    for (noise in noises) {
      for (n in n.obs) {

        print(paste('P:', p, 'N:', n))
        x.params$p <- p
        y.params$noise <- noise
        y.params$thresh <- thresh
        
                                  
        if (y.params$rule == 'mixxor') {
          rr <- paste0(y.params$rule, thresh, 'p_', m.str)
          y.params$thresh1 <- y.params$thresh
          y.params$active1 <- active1
          y.params$k <- length(active1)
          y.params$prob <- mixture.prop

          y.params$thresh2 <- -0.5
          y.params$active2 <- active2
        } else {
          y.params$thresh <- thresh
          rr <- paste0(r, '_', y.params$thresh, '_Identity')
          y.params$active.set <- active1
          y.params$k <- length(active1)
        }

        dir <- paste0('./booleanSims/noise_', noise, '/', rr, '/', 'nObs_', n, '/')
        dir.create(dir, recursive=TRUE)
        iters <- 1:20
        lapply(iters, function(ii) runSim(n, p, x.params, y.params, dir, ii, n.cores))
      }
    }
  }
}



