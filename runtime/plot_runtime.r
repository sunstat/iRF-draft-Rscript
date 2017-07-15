# run runtime_rulefit.r and runtime_irf2.r before running this script

DEBUG = TRUE
pseq = seq(10,80, 10)

runtime_rulefit = rep(0, length(pseq))
auroc_rulefit = rep(0, length(pseq))
for (iter in 1:10){
load(paste0('runtime_rulefit_taly_', iter, '.RData'))
runtime_rulefit = runtime_rulefit+runtime/60
auroc_rulefit = auroc_rulefit + auroc*100
}
runtime_rulefit = runtime_rulefit/10
auroc_rulefit = auroc_rulefit/10

runtime_irf = rep(0, length(pseq))
auroc_irf = rep(0, length(pseq))
for (iter in 1:10){
load(paste0('runtime_irf2_taly_', iter, '.Rdata'))
runtime_irf = runtime_irf+runtime/60
auroc_irf = auroc_irf+auroc[,3]*100
}
runtime_irf = runtime_irf/10
auroc_irf = auroc_irf/10

pseq = seq(10,80, 10)

pdf('plot_runtime_taly_balanced_B_10.pdf', height = 6, width = 12)
par(mfrow=c(1,2))
plot(range(pseq), range(c(runtime_irf, runtime_rulefit))
   , type='n', xlab = 'No. of Features (p)', ylab = 'Runtime (mins.)'
   , main = 'Runtime (Enhancer Data)'
    )
lines(pseq, runtime_irf, type='b', lwd=2, pch=18)
lines(pseq, runtime_rulefit, type='b', lwd=2, pch=18, col='blue')
legend('topleft', legend=c('iRF', 'Rulefit')
     , col=c('black', 'blue'), lwd=2, pch=18, bty='n')

plot(range(pseq), c(75, 85) # range(c(auroc_irf, auroc_rulefit))
   , type='n', xlab = 'No. of Features (p)', ylab = 'AUROC (%)'
   , main = 'AUC (Enhancer Data)'
    )
lines(pseq, auroc_irf, type='b', lwd=2, pch=18)
lines(pseq, auroc_rulefit, type='b', lwd=2, pch=18, col='blue')
legend('topleft', legend=c('iRF', 'Rulefit')
     , col=c('black', 'blue'), lwd=2, pch=18, bty='n')
dev.off()



