DEBUG = TRUE
platform="linux"
rfhome="/persistent/sumbose/iRF_rulefit/rulefit"
source('/persistent/sumbose/iRF_rulefit/rulefit/rulefit.r')
if (!file.exists("/persistent/sumbose/iRF_rulefit/rulefit/akima")){install.packages("akima", lib=rfhome)}
library(akima, lib.loc=rfhome)
require(AUC)

myfun_rulefit <- function(
  x, y
  , xtest, ytest
  , find_interaction = TRUE
  , n_bootstrap = 30
  , cutoff_imp_zscore = 1
  , ... # additional arguments to be passed to rulefit
  ){

  rulefitmode = ifelse(is.factor(y), "class", "regress")

  if (rulefitmode == "class"){
    if ((length(levels(y))!=2))
       stop('only provide binary variables for classification')
    y = as.numeric(levels(y))[y]
    if ((min(y)!=0)|(max(y)!=1))
       stop('only provide 0/1 valued response Y')
    y[y==0]=-1
  }

  rfmod = rulefit(x=x,y=y, rfmode=rulefitmode, ... )

  # predict on test set
  if (!is.null(xtest)){
  yp = rfpred(xtest)
  yp = 1.0/(1.0+exp(-yp))

  print("test set AUC:")  
  print(round(auc(roc(yp, ytest)),2))
  }

  p=ncol(x)
  
  # compute bootstrapped null interaction model 
  cat('compute null models by bootstrap ... ')
  null.models = intnull(ntimes = n_bootstrap, quiet=F)
  cat('done \n')
  
  # find strongly interacting variables
  cat('find strong interacting features ... ')
  int1 = interact(1:p, null.mods=null.models, plot=F)
  imp1var = (int1$int-int1$nullave)/int1$nullstd
  names(imp1var) = 1:p
  select_impvar = which(imp1var > cutoff_imp_zscore)
  p_imp = length(select_impvar)
  store_interaction_1 = imp1var
  cat('done \n')
  print(round(store_interaction_1[select_impvar], 3))
  store_interaction_2 <- NULL
  store_interaction_3 <- NULL
  
  if (p_imp > 1){
  # find interaction strengths of every pair of important variables
  cat('find two order interactions ... ')
  imp2var = array(0, c(p_imp, p_imp))
  rownames(imp2var) = select_impvar
  colnames(imp2var) = select_impvar
  store_interaction_2 = numeric()

  for (i in 1:(p_imp-1)){
     int2 = twovarint(tvar=select_impvar[i]
                    , vars=select_impvar[(i+1):p_imp]
                    , null.mods=null.models, plot=F
                     )
     print(select_impvar[i])

     int_zscore = (int2$int-int2$nullave)/int2$nullstd
     imp2var[i,(i+1):p_imp]=int_zscore
     
     tmp = int_zscore
     names(tmp) = paste(select_impvar[i], select_impvar[(i+1):p_imp], sep="_")
     store_interaction_2 = c(store_interaction_2, tmp)
  }
  
  select_impvar_2 = which(imp2var > cutoff_imp_zscore, arr.ind=TRUE)
  select_impvar_2[,1] = select_impvar[select_impvar_2[,1]]
  select_impvar_2[,2] = select_impvar[select_impvar_2[,2]]
  select_impvar_2_vec = unique(as.vector(select_impvar_2))
  p_imp_2 = length(select_impvar_2_vec)
  print(select_impvar)
  cat('done \n')

  # find strengths of 3 order interactions among every 
  # strongly interacting pairs of variables
  if (p_imp_2 > 2){
  cat('find intereactions of order 3 ... ')
  store_interaction_3 = numeric()
    imp3var = array(0, c(nrow(select_impvar_2), p_imp_2))
    rownames(imp3var) = paste(select_impvar_2[,1], select_impvar_2[,2], sep='_')
    colnames(imp3var) = select_impvar_2_vec

    for (i in 1:nrow(imp3var)){
      target1 = select_impvar_2[i,1]
      target2 = select_impvar_2[i,2]
      others = setdiff(select_impvar_2_vec, c(target1, target2))

      print(c(target1, target2))
      int3 = threevarint(tvar1=target1, tvar2=target2, vars=others
                      , null.mods=null.models, plot=F
                       )
      tmp = (int3$int-int3$nullave)/int3$nullstd
      names(tmp) = paste(paste(target1, target2, sep='_')
                       , others, sep='_'
                        )
      store_interaction_3 = c(store_interaction_3, tmp)

      for (j in 1:length(others))
        imp3var[i, which(colnames(imp3var) == as.character(others[j]))] = tmp[j]

    } # for (i in 1:nrow(imp3var))
  
  cat('done \n')
  } # end if (p_imp > 2)

  } # end if (p_imp > 1)

 out=list()  
 out$yhat = yp
 out$imp1var = imp1var
 store_interaction = store_interaction_1
 if (!is.null(store_interaction_2)){
    store_interaction = c(store_interaction, store_interaction_2)
    out$imp2var = imp2var
 }
 if (!is.null(store_interaction_3)){
    store_interaction = c(store_interaction, store_interaction_3)
    out$imp3var = imp3var
 }
out$interaction = store_interaction

return(out)
}
