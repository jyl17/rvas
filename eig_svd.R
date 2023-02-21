eig_svd <- function(common, k){
  ## return the kth largest svd of the common variants
  ## Args: 
  ## common: absolute path to the common variant matrix
  ## k: number of largest svd one wants to retreive
  ## Output:
  ## svd matrix 
  common <- readRDS(common)
  common <- common[['g']]
  common[is.na(common)] <- 0
  common <- scale(common)
  eig <- svds(common, k)
  u <- eig[['u']]
  u
}