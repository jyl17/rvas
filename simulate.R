power.list=list()

library(MASS)
library(ggplot2)
library(WriteXLS)
install.packages("remotes")
remotes::install_github("grayclhn/dbframe-R-library")
library(dbframe)
library(glmmTMB)
library(ggplot2)
library(mvnfast)
library(Matrix)
# remember dispersion is BCV^2
disps=c(.2,.4,.6,.8,1,2)
#fold changes (added in depletion effects)
fcs=c(1/rev(c(1.01,1.05,1.1,1.2,1.5,2,4)),1.01, 1.05, 1.1, 1.2, 1.5,2,4)
# number of barcode
nreps=1000
# average minimum number of reads per barcode
min.counts=c(100,200,500)
# barcode effect
barcode.variance = c(1,3,5)


eg=expand.grid(disps,fcs,nreps,min.counts, barcode.variance)
colnames(eg)=c('disp', 'fc', 'nbarcode', 'min.count', 'barcode.variance')
eg=data.frame(eg)
eg$estdisp = 0
eg$stddev = 0

start = Sys.time()
for(i in 1:nrow(eg)) {
     n=eg[i,'nbarcode']
     disp=eg[i,'disp']
     min.counts=eg[i,'min.count']
     fold=eg[i,'fc']     
     f=as.factor(c(rep('A',n),rep('B',n)))
     barcode=rbind(diag(1,n), diag(1,n))
     barcode.variance=eg[i,'barcode.variance']
     cov=barcode.variance*(tcrossprod(barcode, barcode))

     # psim=RepParallel(1000, {
     #  # to deal with some stupid random error 
     #  tryCatch({
     #    # annoyingly rnegbin and glm define theta differently from edgeR ... as 1/disp
     #    # use mu and some covariance matrix (symmetric, 1 to 10, based on number of barcodes), 
     #    # to generate multivariate normal and draw new mu from the distribution
     #    # see how that affect theta
     #    # y=MASS::rnegbin(2*n, mu=c(rep(1,n), rep(fold,n))*min.counts, theta=1/disp)
     #    y=MASS::rnegbin(2*n,
     #                    mu=MASS::mvrnorm(mu = c(rep(1,n), rep(fold,n))*min.counts, Sigma = cov),
     #                    theta=1/disp)
     # 
     #    ## assume pre-estimated dispersion, and estimated without error
     #     # gmodel1=glm(y~f, family=MASS::negative.binomial(theta=1/disp)) 
     #     barcodefactor = as.factor(c(rep(1:n),rep(1:n)))
     #     gmodel1 = glmmTMB(y~f + (1|barcodefactor), family = poisson)
     #     p=drop1(gmodel1, test='LRT')[[5]][2]
     #    return(p) 
     #    }, error=function(e) return(NA) )
     # }, simplify='vector', mc.cores=35)
     # psim=na.omit(psim)
     # # note, this is alpha threshold is just a guess
     # eg$power[i]=sum(psim<(.05/5000))/length(psim)
     # print(eg[i,])
     barcode.factor=as.factor(rep(1:n, 2))
     
     # mu1 = exp(MASS::mvrnorm(mu = c(rep(1,n), rep(fold,n)), Sigma = cov))
     diag(cov) <- diag(cov) + 10^(-10)
     mu1 = exp(mvnfast::rmvn(n = 1,mu = c(rep(1,n), rep(fold,n)),sigma = cov))
     y=MASS::rnegbin(2*n,
                     mu=mu1*min.counts,
                     theta=1/disp)
     m1 = glmmTMB(y~f+(1|barcode.factor),family=nbinom2)
     eg$estdisp[i] = 1/sigma(m1)
     temp <- VarCorr(m1)
     eg$stddev[i] = attr(temp$cond$barcode.factor, "stddev")
     print(eg[i,])
}
print(Sys.time() - start)
eg$BCV=as.factor(round(sqrt(eg$disp),2))
eg$disp=as.factor(eg$disp) #round(sqrt(eg$disp),2))

library(ggplot2)
library(viridis)
ggplot(eg, aes(x=fc, y=power, color=disp))+geom_point()+
    geom_line()+
    facet_grid(rows=vars(min.count), cols=vars(nbarcode))+
   # round(,2)
   # [c(1:5,10:14)]
    xlab("fold change")+ scale_x_continuous(trans='log2', breaks=round(fcs,2) , limits=c(.8, 1.2))  +
    #,2), limits=c(.8,1.2)) +
   #  guide = guide_axis(n.dodge=4))+
    scale_color_viridis(discrete=T) +
     #scale_color_brewer(palette="Set1")+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ggtitle('CRISPR editing design power analysis')

x=diag(1,100)
X=rbind(x,x)
XtX=X%*%t(X)
eigen(XtX)$values>0

ggplot(eg, aes(x=disp, y=estdisp)) +
  geom_point(aes(col=factor(min.count))) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(eg, aes(x = barcode.variance, y = stddev^2)) + 
  geom_point(aes(col = factor(fc))) + 
  geom_abline(slope = 1, intercept = 0)
