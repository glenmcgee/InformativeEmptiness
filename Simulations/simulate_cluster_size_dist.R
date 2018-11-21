######################################################
##      Simulate Distribution of Cluster Sizes      ##
######################################################
## Updated: 11/21/2018
## generate dataset with k=1million for each setting and report distribution of cluster sizes
## generating data:
  ## under POISSON + ZIP model
  ## Random intercepts + Random slopes
  ## conditionally + marginally
  ## with and without informativeness

runLOCAL=TRUE ## set to TRUE to run locally ## FALSE is on cluster
set.seed(1000)

#####################
## setting parameters
KK <- 1000000 ## No. of clusters
Xprev <- 0.25 ## Exposure prevalence
a0 <- 0.6; a1 <- -0.2 ## volume model coefficients
b0 <- -4; b1 <- 0.5; b2 <- 0.2 ## outcome model coefficients
e0 <- -2; e1 <- 1 ## zero-inflation model coefficients
#gamma <- -0.1 ## scaling factor for random intercept in volume model
sigma0 <- 2 ## SD of rand eff for X=0
sigma1 <- 1.5 ## SD of rand eff for X=1

###############################################
## set folder to save data and load functions
if(runLOCAL==TRUE){
  path <- "./Simulations/Results/" ## path for results
  funpath <- "./Functions/" ## path for functions
} else{
  path <- "./Simulations/Results/" ## path for results
  funpath <- "./Functions/" ## path for functions
}


####################
## load packages 
source(paste0(funpath,"JMMICS.R"))
source(paste0(funpath,"genICS.R"))
require(glmmML)
require(geepack)
require(JMMICSpack)
require(xtable)



##############################################################
## function to generate conditional and marginal distributions
get_sizeDist <- function(gam=0,m=0,ZIP=FALSE,slopes=FALSE){
  cdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=slopes,cond=TRUE,ZIP=ZIP,   K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma0,sig1.bk=sigma1) ## conditional
  mdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=slopes,cond=FALSE,ZIP=ZIP,  K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma0,sig1.bk=sigma1) ## marginal
  
  cdist <- table(cdat$data_clstlvl$Nk) ## get no. of clusters
  mdist <- table(mdat$data_clstlvl$Nk)
  
  if(m!=0){ ## in case 0s not permitted

    cdist <- c(0,cdist[1:4],sum(cdist[5:length(cdist)])) ## collapse 5+
    mdist <- c(0,mdist[1:4],sum(mdist[5:length(cdist)]))
    
    names(cdist) <- names(mdist) <- c("0","1","2","3","4","5+")  ## rename 
    
  }else{
    
    cdist <- c(cdist[1:5],sum(cdist[6:length(cdist)])) ## collapse 5+
    mdist <- c(mdist[1:5],sum(mdist[6:length(cdist)]))
    
    names(cdist) <- names(mdist) <- c("0","1","2","3","4","5+")  ## rename 
    
  }
  
  cdist <- round(100*cdist/KK) ## get percentages
  mdist <- round(100*mdist/KK)
  
  dists <- data.frame(cbind(mdist,cdist)) ## combine
  
  return(dists)
  
}



##############################################################
## generate distributions for each simulation setting

## Poisson intercepts
dist_pois_gam00_m1 <- get_sizeDist(gam=0.00,m=1,ZIP=FALSE,slopes=FALSE)
dist_pois_gam10_m1 <- get_sizeDist(gam=-0.10,m=1,ZIP=FALSE,slopes=FALSE)
dist_pois_gam25_m1 <- get_sizeDist(gam=-0.25,m=1,ZIP=FALSE,slopes=FALSE)
dist_pois_gam00_m0 <- get_sizeDist(gam=0.00,m=0,ZIP=FALSE,slopes=FALSE)
dist_pois_gam10_m0 <- get_sizeDist(gam=-0.10,m=0,ZIP=FALSE,slopes=FALSE)
dist_pois_gam25_m0 <- get_sizeDist(gam=-0.25,m=0,ZIP=FALSE,slopes=FALSE)

## Poisson slopes
dist_pois_slopes_gam00_m1 <- get_sizeDist(gam=0.00,m=1,ZIP=FALSE,slopes=TRUE)
dist_pois_slopes_gam10_m1 <- get_sizeDist(gam=-0.10,m=1,ZIP=FALSE,slopes=TRUE)
dist_pois_slopes_gam25_m1 <- get_sizeDist(gam=-0.25,m=1,ZIP=FALSE,slopes=TRUE)
dist_pois_slopes_gam00_m0 <- get_sizeDist(gam=0.00,m=0,ZIP=FALSE,slopes=TRUE)
dist_pois_slopes_gam10_m0 <- get_sizeDist(gam=-0.10,m=0,ZIP=FALSE,slopes=TRUE)
dist_pois_slopes_gam25_m0 <- get_sizeDist(gam=-0.25,m=0,ZIP=FALSE,slopes=TRUE)

## ZIP intercepts
dist_ZIP_gam00_m0 <- get_sizeDist(gam=0.00,m=0,ZIP=TRUE,slopes=FALSE)
dist_ZIP_gam10_m0 <- get_sizeDist(gam=-0.10,m=0,ZIP=TRUE,slopes=FALSE)
dist_ZIP_gam25_m0 <- get_sizeDist(gam=-0.25,m=0,ZIP=TRUE,slopes=FALSE)

## ZIP slopes
dist_ZIP_slopes_gam00_m0 <- get_sizeDist(gam=0.00,m=0,ZIP=TRUE,slopes=TRUE)
dist_ZIP_slopes_gam10_m0 <- get_sizeDist(gam=-0.10,m=0,ZIP=TRUE,slopes=TRUE)
dist_ZIP_slopes_gam25_m0 <- get_sizeDist(gam=-0.25,m=0,ZIP=TRUE,slopes=TRUE)



##############################################################
## Make distribution tables

## marginal 
mdist_pois_m1 <- cbind(dist_pois_gam00_m1$mdist,dist_pois_gam10_m1$mdist,dist_pois_gam25_m1$mdist)
mdist_pois_m0 <- cbind(dist_pois_gam00_m0$mdist,dist_pois_gam10_m0$mdist,dist_pois_gam25_m0$mdist)
mdist_pois_slopes_m1 <- cbind(dist_pois_slopes_gam00_m1$mdist,dist_pois_slopes_gam10_m1$mdist,dist_pois_slopes_gam25_m1$mdist)
mdist_pois_slopes_m0 <- cbind(dist_pois_slopes_gam00_m0$mdist,dist_pois_slopes_gam10_m0$mdist,dist_pois_slopes_gam25_m0$mdist)
mdist_ZIP_m0 <- cbind(dist_ZIP_gam00_m0$mdist,dist_ZIP_gam10_m0$mdist,dist_ZIP_gam25_m0$mdist)
mdist_ZIP_slopes_m0 <- cbind(dist_ZIP_slopes_gam00_m0$mdist,dist_ZIP_slopes_gam10_m0$mdist,dist_ZIP_slopes_gam25_m0$mdist)

## conditional 
cdist_pois_m1 <- cbind(dist_pois_gam00_m1$cdist,dist_pois_gam10_m1$cdist,dist_pois_gam25_m1$cdist)
cdist_pois_m0 <- cbind(dist_pois_gam00_m0$cdist,dist_pois_gam10_m0$cdist,dist_pois_gam25_m0$cdist)
cdist_pois_slopes_m1 <- cbind(dist_pois_slopes_gam00_m1$cdist,dist_pois_slopes_gam10_m1$cdist,dist_pois_slopes_gam25_m1$cdist)
cdist_pois_slopes_m0 <- cbind(dist_pois_slopes_gam00_m0$cdist,dist_pois_slopes_gam10_m0$cdist,dist_pois_slopes_gam25_m0$cdist)
cdist_ZIP_m0 <- cbind(dist_ZIP_gam00_m0$cdist,dist_ZIP_gam10_m0$cdist,dist_ZIP_gam25_m0$cdist)
cdist_ZIP_slopes_m0 <- cbind(dist_ZIP_slopes_gam00_m0$cdist,dist_ZIP_slopes_gam10_m0$cdist,dist_ZIP_slopes_gam25_m0$cdist)


mtab <- rbind(  cbind(t(mdist_pois_m1),t(mdist_pois_m0)),
                cbind(t(mdist_pois_slopes_m1),t(mdist_pois_slopes_m0)),
                cbind(matrix(NA,nrow=3,ncol=6),t(mdist_ZIP_m0)),
                cbind(matrix(NA,nrow=3,ncol=6),t(mdist_ZIP_slopes_m0)) )

ctab <- rbind(  cbind(t(cdist_pois_m1),t(cdist_pois_m0)),
                cbind(t(cdist_pois_slopes_m1),t(cdist_pois_slopes_m0)),
                cbind(matrix(NA,nrow=3,ncol=6),t(cdist_ZIP_m0)),
                cbind(matrix(NA,nrow=3,ncol=6),t(cdist_ZIP_slopes_m0)) )

mtab <- cbind(rep(c(0,-0.1,-0.25),times=4),mtab)
ctab <- cbind(rep(c(0,-0.1,-0.25),times=4),ctab)

colnames(mtab) <- colnames(ctab) <- c("Gamma",rep(rownames(dist_pois_gam00_m1),times=2))
rownames(mtab) <- rownames(ctab) <- rep(c("Pois","Pois Slopes","ZIP","ZIP Slopes"),each=3)

xtable(mtab,digits = 0)
xtable(ctab,digits = 0)









##############################################################################
## function to generate conditional and marginal distributions of outcome rate
get_outDist <- function(gam=0,m=0,ZIP=FALSE,slopes=FALSE,dig=1){
  cdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=slopes,cond=TRUE,ZIP=ZIP,   K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma0,sig1.bk=sigma1) ## conditional
  mdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=slopes,cond=FALSE,ZIP=ZIP,  K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma0,sig1.bk=sigma1) ## marginal
  
  cprev <- round(100*mean(cdat$data$Yki),digits=dig) ## get mean outcome
  mprev <- round(100*mean(mdat$data$Yki),digits=dig)
  c0prev <- round(100*mean(cdat$data$Yki[cdat$data$X1ki==0]),digits=dig) ## get mean outcome among unexposed
  m0prev <- round(100*mean(mdat$data$Yki[mdat$data$X1ki==0]),digits=dig)
  c1prev <- round(100*mean(cdat$data$Yki[cdat$data$X1ki==1]),digits=dig) ## get mean outcome among exposed
  m1prev <- round(100*mean(mdat$data$Yki[mdat$data$X1ki==1]),digits=dig)
  


  cN1prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==1]),digits=dig) ## get mean outcome among exposed
  mN1prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==1]),digits=dig)
  cN2prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==2]),digits=dig) ## get mean outcome among exposed
  mN2prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==2]),digits=dig)
  cN3prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==3]),digits=dig) ## get mean outcome among exposed
  mN3prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==3]),digits=dig)
  cN4prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==4]),digits=dig) ## get mean outcome among exposed
  mN4prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==4]),digits=dig)
  cN5prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki>=5]),digits=dig) ## get mean outcome among exposed
  mN5prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki>=5]),digits=dig)
  
  return(list(cprev=cprev,mprev=mprev,
              c01prev=c(c0prev,c1prev),m01prev=c(m0prev,m1prev),
              cNprev=c(cN1prev,cN2prev,cN3prev,cN4prev,cN5prev),
              mNprev=c(mN1prev,mN2prev,mN3prev,mN4prev,mN5prev)))
  
}

##############################################################
## generate outcome prevalences for each simulation setting
## Poisson intercepts
prev_pois_gam00_m1 <- get_outDist(gam=0.00,m=1,ZIP=FALSE,slopes=FALSE)
prev_pois_gam10_m1 <- get_outDist(gam=-0.10,m=1,ZIP=FALSE,slopes=FALSE)
prev_pois_gam25_m1 <- get_outDist(gam=-0.25,m=1,ZIP=FALSE,slopes=FALSE)
prev_pois_gam00_m0 <- get_outDist(gam=0.00,m=0,ZIP=FALSE,slopes=FALSE)
prev_pois_gam10_m0 <- get_outDist(gam=-0.10,m=0,ZIP=FALSE,slopes=FALSE)
prev_pois_gam25_m0 <- get_outDist(gam=-0.25,m=0,ZIP=FALSE,slopes=FALSE)

## Poisson slopes
prev_pois_slopes_gam00_m1 <- get_outDist(gam=0.00,m=1,ZIP=FALSE,slopes=TRUE)
prev_pois_slopes_gam10_m1 <- get_outDist(gam=-0.10,m=1,ZIP=FALSE,slopes=TRUE)
prev_pois_slopes_gam25_m1 <- get_outDist(gam=-0.25,m=1,ZIP=FALSE,slopes=TRUE)
prev_pois_slopes_gam00_m0 <- get_outDist(gam=0.00,m=0,ZIP=FALSE,slopes=TRUE)
prev_pois_slopes_gam10_m0 <- get_outDist(gam=-0.10,m=0,ZIP=FALSE,slopes=TRUE)
prev_pois_slopes_gam25_m0 <- get_outDist(gam=-0.25,m=0,ZIP=FALSE,slopes=TRUE)

## ZIP intercepts
prev_ZIP_gam00_m0 <- get_outDist(gam=0.00,m=0,ZIP=TRUE,slopes=FALSE)
prev_ZIP_gam10_m0 <- get_outDist(gam=-0.10,m=0,ZIP=TRUE,slopes=FALSE)
prev_ZIP_gam25_m0 <- get_outDist(gam=-0.25,m=0,ZIP=TRUE,slopes=FALSE)

## ZIP slopes
prev_ZIP_slopes_gam00_m0 <- get_outDist(gam=0.00,m=0,ZIP=TRUE,slopes=TRUE)
prev_ZIP_slopes_gam10_m0 <- get_outDist(gam=-0.10,m=0,ZIP=TRUE,slopes=TRUE)
prev_ZIP_slopes_gam25_m0 <- get_outDist(gam=-0.25,m=0,ZIP=TRUE,slopes=TRUE)


##############################################################
## Make prevalence tables

## marginal 
mprev_pois_m1 <- c(prev_pois_gam00_m1$mprev,prev_pois_gam10_m1$mprev,prev_pois_gam25_m1$mprev)
mprev_pois_m0 <- c(prev_pois_gam00_m0$mprev,prev_pois_gam10_m0$mprev,prev_pois_gam25_m0$mprev)
mprev_pois_slopes_m1 <- c(prev_pois_slopes_gam00_m1$mprev,prev_pois_slopes_gam10_m1$mprev,prev_pois_slopes_gam25_m1$mprev)
mprev_pois_slopes_m0 <- c(prev_pois_slopes_gam00_m0$mprev,prev_pois_slopes_gam10_m0$mprev,prev_pois_slopes_gam25_m0$mprev)
mprev_ZIP_m0 <- c(prev_ZIP_gam00_m0$mprev,prev_ZIP_gam10_m0$mprev,prev_ZIP_gam25_m0$mprev)
mprev_ZIP_slopes_m0 <- c(prev_ZIP_slopes_gam00_m0$mprev,prev_ZIP_slopes_gam10_m0$mprev,prev_ZIP_slopes_gam25_m0$mprev)

## conditional 
cprev_pois_m1 <- c(prev_pois_gam00_m1$cprev,prev_pois_gam10_m1$cprev,prev_pois_gam25_m1$cprev)
cprev_pois_m0 <- c(prev_pois_gam00_m0$cprev,prev_pois_gam10_m0$cprev,prev_pois_gam25_m0$cprev)
cprev_pois_slopes_m1 <- c(prev_pois_slopes_gam00_m1$cprev,prev_pois_slopes_gam10_m1$cprev,prev_pois_slopes_gam25_m1$cprev)
cprev_pois_slopes_m0 <- c(prev_pois_slopes_gam00_m0$cprev,prev_pois_slopes_gam10_m0$cprev,prev_pois_slopes_gam25_m0$cprev)
cprev_ZIP_m0 <- c(prev_ZIP_gam00_m0$cprev,prev_ZIP_gam10_m0$cprev,prev_ZIP_gam25_m0$cprev)
cprev_ZIP_slopes_m0 <- c(prev_ZIP_slopes_gam00_m0$cprev,prev_ZIP_slopes_gam10_m0$cprev,prev_ZIP_slopes_gam25_m0$cprev)


mprev_tab <- rbind(  cbind((mprev_pois_m1),(mprev_pois_m0)),
                      cbind((mprev_pois_slopes_m1),(mprev_pois_slopes_m0)),
                      cbind(matrix(NA,nrow=3,ncol=1),(mprev_ZIP_m0)),
                      cbind(matrix(NA,nrow=3,ncol=1),(mprev_ZIP_slopes_m0)) )

cprev_tab <- rbind(  cbind((cprev_pois_m1),(cprev_pois_m0)),
                      cbind((cprev_pois_slopes_m1),(cprev_pois_slopes_m0)),
                      cbind(matrix(NA,nrow=3,ncol=1),(cprev_ZIP_m0)),
                      cbind(matrix(NA,nrow=3,ncol=1),(cprev_ZIP_slopes_m0)) )

mprev_tab <- cbind(rep(c(0,-0.1,-0.25),times=4),mprev_tab)
cprev_tab <- cbind(rep(c(0,-0.1,-0.25),times=4),cprev_tab)

colnames(mprev_tab) <- colnames(cprev_tab) <- c("Gamma",rep(c("Prevalence"),times=2))
rownames(mprev_tab) <- rownames(cprev_tab) <- rep(c("Pois","Pois Slopes","ZIP","ZIP Slopes"),each=3)

xtable(mprev_tab,digits = 1)
xtable(cprev_tab,digits = 1)

##############################################################
## Make prevalence tables for exposed and unexposed


## marginal 
m01prev_pois_m1 <- rbind(prev_pois_gam00_m1$m01prev,prev_pois_gam10_m1$m01prev,prev_pois_gam25_m1$m01prev)
m01prev_pois_m0 <- rbind(prev_pois_gam00_m0$m01prev,prev_pois_gam10_m0$m01prev,prev_pois_gam25_m0$m01prev)
m01prev_pois_slopes_m1 <- rbind(prev_pois_slopes_gam00_m1$m01prev,prev_pois_slopes_gam10_m1$m01prev,prev_pois_slopes_gam25_m1$m01prev)
m01prev_pois_slopes_m0 <- rbind(prev_pois_slopes_gam00_m0$m01prev,prev_pois_slopes_gam10_m0$m01prev,prev_pois_slopes_gam25_m0$m01prev)
m01prev_ZIP_m0 <- rbind(prev_ZIP_gam00_m0$m01prev,prev_ZIP_gam10_m0$m01prev,prev_ZIP_gam25_m0$m01prev)
m01prev_ZIP_slopes_m0 <- rbind(prev_ZIP_slopes_gam00_m0$m01prev,prev_ZIP_slopes_gam10_m0$m01prev,prev_ZIP_slopes_gam25_m0$m01prev)

## conditional 
c01prev_pois_m1 <- rbind(prev_pois_gam00_m1$c01prev,prev_pois_gam10_m1$c01prev,prev_pois_gam25_m1$c01prev)
c01prev_pois_m0 <- rbind(prev_pois_gam00_m0$c01prev,prev_pois_gam10_m0$c01prev,prev_pois_gam25_m0$c01prev)
c01prev_pois_slopes_m1 <- rbind(prev_pois_slopes_gam00_m1$c01prev,prev_pois_slopes_gam10_m1$c01prev,prev_pois_slopes_gam25_m1$c01prev)
c01prev_pois_slopes_m0 <- rbind(prev_pois_slopes_gam00_m0$c01prev,prev_pois_slopes_gam10_m0$c01prev,prev_pois_slopes_gam25_m0$c01prev)
c01prev_ZIP_m0 <- rbind(prev_ZIP_gam00_m0$c01prev,prev_ZIP_gam10_m0$c01prev,prev_ZIP_gam25_m0$c01prev)
c01prev_ZIP_slopes_m0 <- rbind(prev_ZIP_slopes_gam00_m0$c01prev,prev_ZIP_slopes_gam10_m0$c01prev,prev_ZIP_slopes_gam25_m0$c01prev)


m01prev_tab <- rbind(  cbind((m01prev_pois_m1),(m01prev_pois_m0)),
                      cbind((m01prev_pois_slopes_m1),(m01prev_pois_slopes_m0)),
                      cbind(matrix(NA,nrow=3,ncol=2),(m01prev_ZIP_m0)),
                      cbind(matrix(NA,nrow=3,ncol=2),(m01prev_ZIP_slopes_m0)) )

c01prev_tab <- rbind(  cbind((c01prev_pois_m1),(c01prev_pois_m0)),
                      cbind((c01prev_pois_slopes_m1),(c01prev_pois_slopes_m0)),
                      cbind(matrix(NA,nrow=3,ncol=2),(c01prev_ZIP_m0)),
                      cbind(matrix(NA,nrow=3,ncol=2),(c01prev_ZIP_slopes_m0)) )

m01prev_tab <- cbind(rep(c(0,-0.1,-0.25),times=4),m01prev_tab)
c01prev_tab <- cbind(rep(c(0,-0.1,-0.25),times=4),c01prev_tab)

colnames(m01prev_tab) <- colnames(c01prev_tab) <- c("Gamma",rep(c("X=0","X=1"),times=2))
rownames(m01prev_tab) <- rownames(c01prev_tab) <- rep(c("Pois","Pois Slopes","ZIP","ZIP Slopes"),each=3)

xtable(m01prev_tab,digits = 1)
xtable(c01prev_tab,digits = 1)



##############################################################
## Make prevalence tables (vs SIZE)

## marginal 
mNprev_pois_m1 <- rbind(prev_pois_gam00_m1$mNprev,prev_pois_gam10_m1$mNprev,prev_pois_gam25_m1$mNprev)
mNprev_pois_m0 <- rbind(prev_pois_gam00_m0$mNprev,prev_pois_gam10_m0$mNprev,prev_pois_gam25_m0$mNprev)
mNprev_pois_slopes_m1 <- rbind(prev_pois_slopes_gam00_m1$mNprev,prev_pois_slopes_gam10_m1$mNprev,prev_pois_slopes_gam25_m1$mNprev)
mNprev_pois_slopes_m0 <- rbind(prev_pois_slopes_gam00_m0$mNprev,prev_pois_slopes_gam10_m0$mNprev,prev_pois_slopes_gam25_m0$mNprev)
mNprev_ZIP_m0 <- rbind(prev_ZIP_gam00_m0$mNprev,prev_ZIP_gam10_m0$mNprev,prev_ZIP_gam25_m0$mNprev)
mNprev_ZIP_slopes_m0 <- rbind(prev_ZIP_slopes_gam00_m0$mNprev,prev_ZIP_slopes_gam10_m0$mNprev,prev_ZIP_slopes_gam25_m0$mNprev)

## conditional 
cNprev_pois_m1 <- rbind(prev_pois_gam00_m1$cNprev,prev_pois_gam10_m1$cNprev,prev_pois_gam25_m1$cNprev)
cNprev_pois_m0 <- rbind(prev_pois_gam00_m0$cNprev,prev_pois_gam10_m0$cNprev,prev_pois_gam25_m0$cNprev)
cNprev_pois_slopes_m1 <- rbind(prev_pois_slopes_gam00_m1$cNprev,prev_pois_slopes_gam10_m1$cNprev,prev_pois_slopes_gam25_m1$cNprev)
cNprev_pois_slopes_m0 <- rbind(prev_pois_slopes_gam00_m0$cNprev,prev_pois_slopes_gam10_m0$cNprev,prev_pois_slopes_gam25_m0$cNprev)
cNprev_ZIP_m0 <- rbind(prev_ZIP_gam00_m0$cNprev,prev_ZIP_gam10_m0$cNprev,prev_ZIP_gam25_m0$cNprev)
cNprev_ZIP_slopes_m0 <- rbind(prev_ZIP_slopes_gam00_m0$cNprev,prev_ZIP_slopes_gam10_m0$cNprev,prev_ZIP_slopes_gam25_m0$cNprev)


mNprev_tab <- rbind(  cbind((mNprev_pois_m1),(mNprev_pois_m0)),
                cbind((mNprev_pois_slopes_m1),(mNprev_pois_slopes_m0)),
                cbind(matrix(NA,nrow=3,ncol=5),(mNprev_ZIP_m0)),
                cbind(matrix(NA,nrow=3,ncol=5),(mNprev_ZIP_slopes_m0)) )

cNprev_tab <- rbind(  cbind((cNprev_pois_m1),(cNprev_pois_m0)),
                cbind((cNprev_pois_slopes_m1),(cNprev_pois_slopes_m0)),
                cbind(matrix(NA,nrow=3,ncol=5),(cNprev_ZIP_m0)),
                cbind(matrix(NA,nrow=3,ncol=5),(cNprev_ZIP_slopes_m0)) )

mNprev_tab <- cbind(rep(c(0,-0.1,-0.25),times=4),mNprev_tab)
cNprev_tab <- cbind(rep(c(0,-0.1,-0.25),times=4),cNprev_tab)

colnames(mNprev_tab) <- colnames(cNprev_tab) <- c("Gamma",rep(c("Nk=1","2","3","4","5+"),times=2))
rownames(mNprev_tab) <- rownames(cNprev_tab) <- rep(c("Pois","Pois Slopes","ZIP","ZIP Slopes"),each=3)

xtable(mNprev_tab,digits = 1)
xtable(cNprev_tab,digits = 1)
