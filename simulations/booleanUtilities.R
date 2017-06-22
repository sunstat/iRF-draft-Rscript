runSim <- function(n, p, x.params, y.params, dir, seed, n.cores=1) {
  # Wrapper function to run simulation with given input parameters.
  # args:
  #   n: number of observations
  #   p: number of features
  #   x.params: named list containing feature generating parameters
  #   y.params: named list containing response generating parameters
  #   dir: output directory where simulations should be saved
  #   seed: seed for random number generation
  #   n.cores: number of cores used in iRF
  set.seed(seed)
  x <- genX(n, x.params)
  y <- genY(x, y.params)
  ns <- ifelse(is.factor(y), 1, 10)
  if (!is.null(x.params$cov)) x.params$cov <- NULL
  
  int.return <- 1:5
  fit <- iRF(y=y, x=x, 
             interactions.return=int.return,
             n.iter=max(int.return),
             nodesize=ns,
             n.bootstrap=20,
             ntree=500,
             n.core=n.cores)

  gini.x <- decreaseGiniWrapper(x, y)
  file <- paste0(dir, 'iter', seed, '.Rdata')
  save(file=file, fit, gini.x, x.params, y.params)

}

genX <- function(n, x.params) {
  # Generates n x p design matrix 
  # args:
  #   n: number of rows (observations) in design matrix
  #   x.params: named list specifying data generating parameters including:
  #     <dist> (distribution of x), <p> (number of features), <struct>
  #     (structure of covariance matrix, optional), <rho> (optional)
  p <- x.params$p
  require(LaplacesDemon)
  
  if (x.params$dist == 'bernoulli')
    x <- matrix(rbinom(n*p, size=1, prob=x.params$prob), nrow=n)
  else if (x.params$dist == 'cauchy' & is.null(x.params$struct))
    x <- matrix(rcauchy(n*p), nrow=n)
  else {
    cov.mat <- makeMat(p, x.params$struct, x.params$rho, k=10)
    x <- rmvc(n=n, mu=rep(0, p), S=cov.mat$S)
  }
  return(x)
}

genY <- function(x, y.params) {
  # Generate response in simulations of iRF paper
  # args:
  #   x: design matrix
  #   y.params: named list specifying response generating parameters parameters
  #     inluding: <rule> (generative model one of or, and, xor, mixxor),
  #     <thresh> (thresold for generative model), <active> (active features),
  #     <classification> (TRUE/FALSE indicating whether classification or
  #     regression iRF), <k> (number of distinct xor components if rule is xor),
  #     <prob> (mixture propoortions if rule is mixxor)
  n <- nrow(x)
  p <- ncol(x)
 
  y <- switch(y.params$rule,
              or=genOR(x, y.params$thresh, y.params$active),
              and=genAND(x, y.params$thresh, y.params$active),
              xor=genXOR(x, y.params$thresh, y.params$active, y.params$k),
              mixxor=genMixXOR(x, y.params$prob, y.params$thresh1, 
                               y.params$thresh2, y.params$active1, 
                               y.params$active2, y.params$k))
         
  if (y.params$noise == 0) {
    y <- as.factor(as.numeric(y))
  } else if (y.params$classification) {
    n.class <- trunc(table(y) * y.params$noise)
    idcs.swap <- mapply(function(cc, nn)
                        sample(which(y == cc), nn, replace=TRUE),
                        as.factor(names(n.class)), n.class)
    idcs.swap <- unlist(idcs.swap)
    y[idcs.swap] <- ifelse(y[idcs.swap], FALSE, TRUE)
    y <- as.factor(as.numeric(y))
  } else {
    y <- as.numeric(y) + rnorm(n, sd=y.params$noise)
  }
  return(y)
}

genOR <- function(x, tt, active.set=1:4) {
  # Generate response under "OR" rule between active variables
  # args:
  #   x: data matrix
  #   tt: threshold, x > tt are "active"
  #   active set: active set over which response will be generated

  stopifnot(ncol(x) >= max(active.set))
  y <- apply(x[,active.set] > tt, MAR=1, any)
  return(y)
}

genAND <- function(x, tt, active.set) {
  # Generate response under "AND" rule between active variables
  # args:
  #   x: data matrix
  #   tt: threshold, x > tt are "active"
  #   active set: active set over which response will be generated
  stopifnot(ncol(x) >= max(active.set))
  y <- apply(x[,active.set] > tt, MAR=1, all)
  return(y)
} 


genXOR <- function(x, tt, active.set, k=NULL) {
  # Generate response under XOR simulation
  # args:
  #   x: data matrix
  #   tt: threshold, x > tt are "active"
  #   active set: active set over which response will be generated
  #   type: what type of xor should be generated
  #   k: number of distinct xor components
  stopifnot(ncol(x) >= max(active.set))

  # standard XOR data generating mechanism
  stopifnot(length(active.set) %% k == 0)
  active.grouped <- matrix(active.set, ncol=k, byrow=TRUE)
  xorGrouped <- function(p) rowSums(x[,p] > tt) %% 2 == 1
      
  group.resp <- apply(active.grouped, MAR=1, xorGrouped)
  y <- apply(as.matrix(group.resp), MAR=1, any)
  return(y)
}


genMixXOR <- function(x, prop, tt1, tt2, active1, active2, k) {
  # Generate responses under mixture of xor/and rules
  # args:
  #   x: data matrix
  #   prop: mixture proportions
  #   tt1: threshold for xor component
  #   tt2: threshold for and component
  #   active1: active features for xor component
  #   active2: active features for and component
  #   k: number of distinct xor components for xor component
   
  x1.idcs <- sample(nrow(x), trunc(p * nrow(x)))
  x1 <- x[x1.idcs,]
  x2 <- x[-x1.idcs,]

  y.xor <- genXOR(x1, tt1, active1, k=k)
  y.and <- genAND(x2, tt2, active2)
  
  y <- rep(FALSE, nrow(x))
  y[x1.idcs] <- y.xor
  y[-x1.idcs] <- y.and
  return(y)
}


gini <- function(y) {
  # Compute the two class gini impurity for responses
  if (is.factor(y)) yy <- as.numeric(y) - 1
  else yy <- y
  return(1 - mean(yy) ^ 2 - mean(1 - yy) ^ 2)
}

decreaseGini <- function(x, y) {
  # Compute the decrease in gini based on best initial split for a single variable x 
  n <- length(y)
  gini.y <- gini(y)
  x.order <- order(x)
  y.sort <- y[x.order]
  split.gini <- sapply(1:(n-1), function(i) 
    c(gini(y.sort[1:i]), gini(y.sort[-(1:i)])))
  weights <- rbind(1:(n-1), (n-1):1)
  child.gini <- colSums(weights * split.gini) / n
  dec.gini <- gini.y - min(child.gini)

  return(c(y.gini=gini.y, child.gini=min(child.gini), dec.gini=dec.gini))
}


decreaseGiniWrapper <- function(x, y) {
  # Compute the decrease in gini over all variables in x
  gg <- apply(x, MAR=2, decreaseGini, y)
  return(gg)
}

