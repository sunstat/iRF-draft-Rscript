library(iRF)
load('../data/enhancer.Reata')
load('../data/rfSampleSplitNoise.Rdata')

tf.idcs <- 46:80
X <- X[,tf.idcs]
varnames.agg <- varnames.all[tf.idcs, 2]
varnames.agg <- gsub('_', '.', varnames.agg)

responseTF <- function(x) {
  # Generate synthetic response from enhancer data 
  y <- (x[,'kr1'] > 1.25 | x[,'kr2'] > 1.25) & 
    (x[,'twi1'] > 1 | x[,'twi2'] > 1) & 
    (x[,'D1'] > 1.25) & 
    (x[,'hb1'] > 1.5 | x[,'hb2'] > 1.5) & 
    (x[,'wt_ZLD'] > 2)
    return(y)
}


n <- nrow(X)
p <- ncol(X)
n.reps <- 20
sample.sizes <- c(200, 400, 600, 800, 1000)
noise <- c('rf', 'swap')
dir <- 'simulations/tf/'
dir.create(dir, recursive=TRUE)

for (ns in noise) {
  for (ss in sample.sizes) {
    write.dir <- paste0(dir, 'noise_', ns, '/', 'nTrain_', ss, '/')
    dir.create(write.dir, recursive=TRUE)

    lapply(1:n.reps, function(i) {

      # Generate responses and introduce additional interactions
      if (ns == 'rf') {
        pred.prob <- pmin(2 * apply(pred.prob, 1, max), 1)
        y.noise <- as.numeric(runif(n) < pred.prob)
        y <- responseTF(x=X)
        
        # take balanced sample of responses
        n1 <- sum(y)
        n0 <- sum(!y)
        i1 <- sample(which(y), min(c(n1, n0)))
        i0 <- sample(which(!y), min(c(n1, n0)))
        y <- y | y.noise
        x.balanced <- X[c(i1, i0),]
        y <- y[c(i1, i0)]

      } else if (ns == 'swap') {

        y <- responseTF(x=X)
        
        # take balanced sample of responses
        n1 <- sum(y)
        n0 <- sum(!y)
        i1 <- sample(which(y), min(c(n1, n0)))
        i0 <- sample(which(!y), min(c(n1, n0)))
        x.balanced <- X[c(i1, i0),]
        y <- y[c(i1, i0)]
        idcs.flip <- sample(length(y), round(length(y) * 0.2))
        if (length(idcs.flip)) y[idcs.flip] <- !y[idcs.flip]
      } 

      idcs.balance <- c(i1, i0)
      prop.1 <- mean(y)
      prop.0 <- 1 - prop.1
      train.id <- c(sample(which(y), round(ss * prop.1), replace=FALSE)
                  , sample(which(!y), round(ss * prop.0), replace=FALSE))
      test.id <- setdiff(1:nrow(x.balanced), train.id) 
      y <- as.factor(as.numeric(y))

      fit <- iRF(x=x.balanced[train.id,], y=y[train.id], 
                 xtest=x.balanced[test.id,], ytest=y[test.id], 
                 n.iter=5, interactions.return=1:5, 
                 n.bootstrap=20, varnames.grp=varnames.agg,
                 n.cores=n.cores) 

      y.test <- y[test.id]
      y.train <- y[train.id]
      train.id <- idcs.balance[train.id]
      test.id <- idcs.balance[test.id]
      save(file=paste0(write.dir, 'rep', i, '.Rdata'), fit, train.id, test.id, y.train, y.test)
    })

  }
}
