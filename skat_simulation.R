skat_simulation <- function(pC, Z, gvec.cor, causal_prop, error_sd, len, genenum){
  ## simulation based on specified proportion or causal variant and random effect
  ## Args: 
  ## pC: absolute path to variant annotation file
  ## Z: absolute path to dosage matrix
  ## gvec.cor: absolute path to correlation matrix
  ## causal_prop: a vector of causal prop (from 0.0 to 1.0) 
  ## error_sd: a vector of random effect (from 0.0 to 2.0)
  ## len: length of the vector of causal_prop
  ## genenum: number of genes of interest
  ## Output:
  ## matrix of p value for each gene with p value
  ## to retreive the power, here's an example with 
  ## calsal_prop = seq(0.1,0.5,length.out = 5) and error_sd = seq(0.3,1.6,length.out = 5)
  ## causal_prop = seq(0.1,0.5,length.out = 5)
  ## causal_prop = rep(causal_prop, each = 5)
  ## causal_prop = factor(causal_prop)
  ## error.sd = seq(0.3,1.6,length.out = 5)
  ## error.sd = rep(error.sd, times = 5)
  ## power = vector(length = 25)
  ## for (n in 1:25){
  ##   temp = pvalue[[n]]
  ##   power[n] = sum(temp < 1e-4)/100
  ## }
  
  
  # extract all the submatrices per gene for variants that don't have synonymous coding consequences 
  nonsyn.matrices=list()
  
  # extract all the variants per chromosome with minor allele > 5%
  common.variants=list()
  
  pC=readRDS(pC)
  Z=readRDS(Z)
  
  #calculate variant allele frequency 
  af=colSums(Z)/(nrow(Z)*2)
  #convert to minor allele frequency, make sure that the least common one at a frequency less than 0.5
  af[af>.5]=(1-af[af>.5])
  
  common.variants[[chr]]=Z[, af>.05]
  
  pC.subset=pC[pC$CONSEQUENCE!='synonymous',]
  
  #split out all the nonsynonymous and nonsense variants 
  spC=split(pC.subset, pC.subset$GENEID)
  
  
  for(gene in names(spC)){
    m = colnames(Z) %in% spC[[gene]]$idx
    # subset information about the nonsynonymous variants
    nonsyn.matrices[[gene]]=Z[,m]
  }
  
  gvec.cor<- readRDS(gvec.cor)
  
  causal_prop = causal_prop
  proportion.causal = vector()
  sd.error = vector()
  power.total = vector()
  pvalue = vector()
  
  for (n in causal_prop){
    temp = rep(n, len)
    proportion.causal = append(proportion.causal, temp)
    error_sd = error_sd
    for (error in error.sd){
      sd.error = append(sd.error, error)
      power = 0
      pv = vector(length = 10)
      for (j in 1:genenum){
        gene=names(nonsyn.matrices)[j]
        Z= nonsyn.matrices[[gene]]
        maf= colSums(Z)/(nrow(Z)*2)
        rar=(which(maf<.01))
        ncaus=round(length(rar)*n)
        B=rep(0, length(maf))
        B[sort(sample(rar, ncaus))]=-1
        #get predicted strain effect
        XB=Z%*%B
        
        #flooring operation for strains with multiple variants
        XB[XB<0]=-1
        
        simy=XB+rnorm(nrow(gvec.cor),mean=0, sd=error)
        X = matrix(1, length(simy))
        obj = SKAT_NULL_emmaX(simy~ X, K = gvec.cor)
        test=SKAT(Z,obj, is_dosage=T, is_check_genotype=F)
        pv[j] = test$p.value
      }
      pv = list(pv)
      pvalue = append(pvalue, pv)
    }
  }
  
}