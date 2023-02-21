skat_joint <- function(pheno, gvec.cor, unique.chrs, u, our.dir){
  ## SKAT without effect of common variant
  ## Args:
  ## pheno: absolute path to phenotype matrix with rowname variant and column name as phenotype
  ## gvec.cor: absolute path to covariance matrix among all strains 
  ## (correlation matrix based on both common and rare variants)
  ## unique.chrs: names of chromosomes one wants to investigate
  ## u: absolute path to the svd (see eig_svd function)
  ## out.dir: absolute path to the directory that contains variant annotation file 
  ## (only coding genetic variants) and dosage matrix
  ## Output: matrix of p value with rowname phenotype and column name gene
  
  
  pheno <- readRDS(pheno)
  gvec.cor<- readRDS(gvec.cor)
  u <- readRDS(u)
  u <- u[which(rownames(gvec.cor) %in% rownames(pheno)),]
  # get rid of extra entries in the correlation matrix
  gvec.cor <- gvec.cor[which(rownames(gvec.cor) %in% rownames(pheno)), ]
  gvec.cor <- gvec.cor[, which(colnames(gvec.cor) %in% rownames(pheno))]
  name <- rownames(pheno)
  out.dir = out.dir
  unique.chrs = unique.chrs
  datalist = vector("list")
  
  for (i in colnames(pheno)){
    print(i)
    y <- pheno[[i]]
    y <- setNames(y, name)
    X <- matrix(1, length(y))
    X <- cbind(X,u)
    obj <-  SKAT_NULL_emmaX(y~ X, K = gvec.cor)
    p <- c()
    genename <- c()
    for (chr in unique.chrs){
      print(chr)
      # extract all the submatrices per gene for variants that don't have synonymous coding consequences 
      nonsyn.matrices=list()
      
      # extract all the variants per chromosome with minor allele > 5%
      common.variants=list()
      rare.variants=list()
      
      pC=readRDS(file=paste0(out.dir , 'pC_', chr, '.RDS'))
      Z=readRDS(file=paste0(out.dir, 'Z_', chr, '.RDS'))
      
      #calculate variant allele frequency 
      af=colSums(Z)/(nrow(Z)*2)
      #convert to minor allele frequency, make sure that the least common one at a frequency less than 0.5
      af[af>.5]=(1-af[af>.5])
      
      pC.subset=pC[pC$CONSEQUENCE!='synonymous',]
      
      #split out all the nonsynonymous and nonsense variants 
      spC=split(pC.subset, pC.subset$GENEID)
      
      
      for(gene in names(spC)){
        m = colnames(Z) %in% spC[[gene]]$idx
        # subset information about the nonsynonymous variants
        nonsyn.matrices[[gene]]=Z[,m]
      }
      
      for (j in 1:length(nonsyn.matrices)){
        gene=names(nonsyn.matrices)[j]
        Z = nonsyn.matrices[[gene]]
        rare <- names(af[colnames(Z)] < 0.05)
        Z <- Z[, rare]
        Z <- Z[which(rownames(Z) %in% rownames(pheno)), ]
        pvalue = tryCatch(
          {
            test <- SKAT(Z, obj, is_dosage = T, is_check_genotype = F)
            test$p.value
          },
          error = function(e){
            NA
          }
        )
        p <- append(p, pvalue)
        genename <- append(genename, gene)
      }
    }
    p <- setNames(p, genename)
    datalist[[i]] <- p
  }
  big_data = do.call(rbind, datalist)
  big_data
}