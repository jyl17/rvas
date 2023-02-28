library(BEDMatrix)
library(rrBLUP)
library(NAM)
library(tidyverse)
library(regress)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(VariantAnnotation)
library(GenomicFeatures)
library(vcfR)
library(SKAT)

# test on real data
u <- eig_svd("~/Desktop/Research/GWAS/rvas/yeast_common_variants.RDS", 10)

data <- skat_joint("~/Desktop/Research/GWAS/rvas/pheno.RDS", 
                   "~/Desktop/Research/GWAS/rvas/data/all_nonsyn_cor.RDS", 
                   paste0('chr', as.roman(1:16)),
                   u, '~/Desktop/Research/GWAS/rvas/ref/')

# simulation
# simulated p value for different combo of sd and random effect
sim <- skat_simulation('~/Desktop/Research/GWAS/rvas/ref/PC_chrI.RDS', 
                '~/Desktop/Research/GWAS/rvas/ref/Z_chrI.RDS',
                "~/Desktop/Research/GWAS/rvas/data/all_nonsyn_cor.RDS",
                seq(0.1,0.5,length.out = 5),
                seq(0.3,1.6,length.out = 5),
                5, 50)

