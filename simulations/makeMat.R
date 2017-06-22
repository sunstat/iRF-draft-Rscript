#################################################################
# Wrapper function for producing many kinds of covariance matrices
##################################################################


makeMat <- function(p, 
                    method = c('identity', 'toeplitz', 'exp.decay', 
                               'equicorr', 'sp.banded', 'sp.random', 
                               'dense', 'img.banded', 'block'),
                    rho = NULL, 
                    sp = NULL, 
                    lowe = NULL, 
                    hie = NULL, 
                    nrow = NULL,
                    k=NULL){
  
  # a function for producing several kinds of covariance matrices
  # toeplitz, exp.decay and equicorr as in Dezeure (hdi) paper
  # requires one parameter, rho (default values as in Dezeure-hdi)
  # sp.banded produces a banded sparse matrix
  # sp.random produces a sparse matrix with randomly placed entries
  # both requires three parameters: sp (sparsity factor), 
  #   lowe (minimum eigenvalue), hie (maximum eigenvalue)
  # dense produces a dense matrix
  # requires lowe (minimum eigenvalue), hie (maximum eigenvalue)
  #
  # args:
  #   p: dimension of covariance matrix
  #   method: what kind of coavariance matrix to construct
  #   rho: correlation
  #   k: number of blocks in block covariance matrix
  
  # output contains four objects:
  # S, the covariance matrix
  # Sinv, its inverse
  # Sroot, its square root
  # elist, its list of eigenvalues
  if (method =='identity') {
    S <- diag(p)
    inverse.return <- FALSE
  }
  
  # toeplitz construction
  if (method == 'toeplitz'){
    if(is.null(rho)){
      rho <- 0.9
    }
    
    init.row <- rho ^ abs(1 - 1:p)
    S <- toeplitz(init.row) 
    inverse.return = F
  } #end of toeplitz
  
  # exp.decay construction
  if (method == 'exp.decay'){
    if(is.null(rho)){
      rho = 0.4
    }
    
    S = matrix(0, p, p)
    for (i in 1:p){
      for (j in i:p){
        temp <- rho ^ (abs(i - j) / 5)
        S[i, j] <- temp
        S[j, i] <- temp
      }
    }
    
    inverse.return = T
  } #end of exp.decay
  
  # equicorr construction
  if(method == 'equicorr'){
    if(is.null(rho)){
      rho = 0.8
    }
    
    S = (1-rho) * diag(p) + rho * matrix(1, p, p)
    
    inverse.return = F
  } # end of equicorr
  
  # sparse banded construction
  if(method == 'sp.banded'){
    if(any(is.null(c(sp, lowe, hie)))){
      print('Parameters not specified')
    } else {
      
      # fixing limits on eigenvalues, because inverse is returned
      lowei = 1/hie
      hiei = 1/lowe
      
      nums = sort(runif(sp,1/(sp+1),1), decreasing = T) * sample(c(-1,1), sp, replace = T)
      S = matrix(0, p, p)
      for (i in 1:(p-1)){
        istart = i+1
        iend = min(i+sp,p)
        ilength = iend - istart + 1
        S[i,istart:iend] = nums[1:ilength]
        S[istart:iend,i] = nums[1:ilength]
      }
      es = eigen(S)$values
      mine = min(es)
      maxe = max(es)
      S = (hiei - lowei)/(maxe - mine) * S + (lowei - (hiei - lowe)*mine/(maxe-mine)) * diag(p)
      
      inverse.return = T
    }
  } # end of sp.banded
  
  # sparse random construction
  if (method == 'sp.random'){
    if(any(is.null(c(sp, lowe, hie)))){
      print('Parameters not specified')
    } else {
      
      # fixing limits on eigenvalues, because inverse is returned
      lowei = 1/hie
      hiei = 1/lowe
      
      S = matrix(0, p, p)
      for (i in 1:p){
        u = array(0, dim = p)
        u[sample(p,sp)] = rnorm(sp)
        S = S + tcrossprod(u)
      }
      es = eigen(S)$values
      mine = min(es)
      maxe = max(es)
      S = (hiei - lowei)/(maxe - mine) * S + (lowei - (hiei - lowei)*mine/(maxe-mine)) * diag(p)
      inverse.return = T
    }
  } # end of sp.random
  
  # dense construction
  if (method == 'dense'){
    if(any(is.null(c(lowe, hie)))){
      print('Parameters not specified')
    } else {
      # fixing limits on eigenvalues, because inverse is returned
      lowei = 1/hie
      hiei = 1/lowe
      
      S = matrix(0,p,p)
      for (i in 1:p){
        uvec = rnorm(p)
        S = S + tcrossprod(uvec)
      }
      es = eigen(S)$values
      mine = min(es)
      maxe = max(es)
      S = (hiei - lowei)/(maxe - mine) * S + (lowei - (hiei - lowei)*mine/(maxe-mine)) * diag(p)
      
      inverse.return = T
    }
  } # end of dense
  
  if (method == 'block') {
    if (is.null(k)) stop('NUMBER OF BLOCKS NOT SPECIFIED')
    S <- blockMat(p, k, rho)
    inverse.return <- FALSE
  }
  
  if(inverse.return){
    Sinv = S
    esi = eigen(Sinv)
    S = esi$vectors %*% diag(1/esi$values) %*% t(esi$vectors)
    Sroot = esi$vectors %*% diag(1/sqrt(esi$values)) %*% t(esi$vectors)
    SinvRoot = esi$vectors %*% diag(sqrt(esi$values)) %*% t(esi$vectors)
    elist = sort(1/esi$values, decreasing = T)
  } else {
    es = eigen(S)
    Sinv = es$vectors %*% diag(1/es$values) %*% t(es$vectors)
    SinvRoot = es$vectors %*% diag(sqrt(1/es$values)) %*% t(es$vectors)
    Sroot = es$vectors %*% diag(sqrt(es$values)) %*% t(es$vectors)
    elist = es$values
  }
  
  
  res <- list()
  res$S <- S
  res$Sinv <- Sinv
  res$SinvRoot <- SinvRoot
  res$Sroot <- Sroot
  res$elist <- elist
  return(res)
}

blockMat <- function(p, k, rho) {
  # p: dimension of matrix
  # k: number of blocks
  # rho: within block correlation, everthing else 0
  stopifnot(p %% k == 0)
  blk.sz <- p / k
  blocks <- lapply(1:k, function(i) (1 + (i - 1) * blk.sz):(i * blk.sz))
  
  sigma <- matrix(0, nrow=p, ncol=p)
  for (b in blocks) sigma[b,b] <- rho
  diag(sigma) <- 1
  return(sigma)
} 
