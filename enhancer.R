library(iRF)
n.cores <- 24
set.seed(47)
load('./data/enhancer.Rdata')

fit <-  iRF(x=X[train.id,], 
           y=as.factor(Y[train.id]), 
           xtest=X[test.id,], 
           ytest=as.factor(Y[test.id]), 
           n.iter=3,
           interactions.return=1:3, 
           varnames.grp=varnames.all$Predictor_collapsed, 
           n.core=n.cores, 
           n.bootstrap=30
          )
save(file='results_enhancer.Rdata', fit)


