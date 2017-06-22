library(iRF)
load('./data/splicing.Rdata')
n.cores <- 24
set.seed(47)

fit <- iRF(x=x[train.id,], y=y[train.id],
           n.iter=3,
           interactions.return=1:3,
           n.bootstrap=1, 
           n.core=n.cores, 
           varnames.grp=varnames.all$Predictor_collapsed)
save(file='results_splice.Rdata', fit)
