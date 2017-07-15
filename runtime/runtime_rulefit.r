# first change directory of rulefit programs in ../rulefit/library_rulefit.r 
rm(list=ls(all=TRUE))
DEBUG = TRUE
source('rulefit/library_rulefit.r')
require(doMC)
registerDoMC()
options(cores=15)
set.seed(32)

load('Y_X_Taly.RData')
drop_vars = c('ZLDm_H3', 'ZLDm_H3K18ac', 'ZLDm_H3K4me1', 'ZLDm_ZLD')
X = X[,!(colnames(X) %in% drop_vars)]
load('varnames_all.RData')
varnames_all = varnames_all[!(varnames_all[,1] %in% drop_vars),]
varnames_agg = varnames_all[,2]

rownames(X) <- NULL
n0=sum(Y==0); n1=sum(Y==1)
if (DEBUG){
train.id=c(sample(which(Y==1), round(n1/2)), sample(which(Y==0), round(n1/2)))
test.id=c(setdiff(which(Y==1), train.id), sample(setdiff(which(Y==0), train.id), round(n1/2)))
}

for(iter in 1:10){
p = ncol(X)
col_permuted= sample(1:p, p, replace=FALSE)
print(col_permuted)
pseq = seq(10, p, by=10)

runtime = rep(0, length(pseq))
auroc = rep(0, length(pseq))

for (i in 1:length(pseq)){
print(paste('p --->', pseq[i]))
id =  col_permuted[1:pseq[i]]
ptm <- proc.time()

ff = myfun_rulefit(x=X[train.id,id]
       , y=as.factor(Y[train.id])
       , xtest=X[test.id,id]
       , ytest=as.factor(Y[test.id])
       , find_interaction = TRUE
       , n_bootstrap = 10
       , cutoff_imp_zscore = 1
         )

runtime[i] = (proc.time()-ptm)[3]
auroc[i] = auc(roc(ff$yhat, as.factor(Y[test.id])))
save(runtime, auroc, file=paste0('runtime_rulefit_taly_', iter, '.RData'))
}
}

rm(list=ls(all=TRUE))



