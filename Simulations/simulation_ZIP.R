##############################
##      ICS Simulation      ##
##############################
## Updated: 11/21/2018
## generating data:
  ## under ZIP model
  ## conditionally and marginally
  ## with and without informativeness

runLOCAL=TRUE ## set to TRUE to run locally ## FALSE is on cluster

#####################
## setting parameters
RR <- 3000 ## No. of iterations
KK <- 4000 ## No. of clusters
JJ <- 300 ## Number of jobs in array (only matters on cluster)
Xprev <- 0.25 ## Exposure prevalence
a0 <- 0.6; a1 <- -0.2 ## volume model coefficients
b0 <- -4; b1 <- 0.5; b2 <- 0.2 ## outcome model coefficients
e0 <- -2; e1 <- 1 ## zero-inflation model coefficients
#gamma <- -0.1 ## scaling factor for random intercept in volume model
sigma <- 2 ## standard deviation of random intercepts
outputsize <- 85 ## width of output (ie number of variables reported) DO NOT CHANGE UNLESS MODELS/FITS ARE CHANGED

###############################################
## set folder to save data and load functions
if(runLOCAL==TRUE){
  path <- "../Simulations/Results/" ## path for results
  funpath <- "../Functions/" ## path for functions
} else{
  path <- "../Simulations/Results/" ## path for results
  funpath <- "../Functions/" ## path for functions
}


####################
## load packages 
source(paste0(funpath,"JMMICS.R"))
source(paste0(funpath,"genICS.R"))
require(glmmML)
require(geepack)
require(JMMICSpack)

####################
## set up parallel
if(runLOCAL==TRUE){
  loopvec <- 1:RR ## run entire loop
  suffix <- "" ## output file doesnt have a suffix
}else{
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  iter_no <- as.integer(args[1])  ## here there is only one argument, the iteration number
  loopvec <- ((RR/JJ)*(iter_no-1) +1):((RR/JJ)*iter_no) ## portion of the loop to run
  suffix <- paste0("_iter",iter_no) ## append output file names with the iteration number, to be combined later
}








#############################
##  Fit conditional models
fit_mods_cond <- function(gendat,size_plus=0,ZIP=TRUE){
  df <- gendat$data; df2 <- gendat$data_clstlvl
  
  ## Naive GLMM
  # (glmer doesnt give SE for sigma)
  naive <- glmmML(Yki~X1ki + X2ki,cluster=IDki,family=binomial,data=df,method="ghq",n.points=60, control=list(maxit=2000)) 
  est_naive <- c(naive$coefficients,naive$sigma); SE_naive <- c(naive$coef.sd,naive$sigma.sd)
  if(length(SE_naive)<length(est_naive)){SE_naive <- rep(NA,length(est_naive))} ## if no convergence
  MSEP_naive <- mean((naive$posterior.modes-gendat$true_randeff[df2$Nk!=0])^2)
  
  ## Joint - No 0s
  ## NOTE: since we are only considering the non-empty clusters, we assume a POISSON model for this analysis
  IDk_non0 <- which(df2$Nk!=0) ## which clusters are nonempty
  jno0s <- JMMICS_fit(Nk=df2$Nk[IDk_non0],Zk=cbind(1,df2$X1k)[IDk_non0,],Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=as.numeric(as.factor(df$IDki)),Z0k=cbind(1,df2$X1k),minNk=size_plus,ZIP=FALSE,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad = 60)
  est_jno0s <- jno0s[[1]]; SE_jno0s <- jno0s[[2]]
  if(length(SE_jno0s)<length(est_jno0s)){SE_jno0s <- rep(NA,length(est_jno0s))} ## if no convergence
  MSEP_jno0s <- mean((jno0s[[3]]-gendat$true_randeff[df2$Nk!=0])^2)
  
  ## Joint Model
  joint <- JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=df$IDki,Z0k=cbind(1,df2$X1k),minNk=size_plus,ZIP=ZIP,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad = 60)
  est_joint <- joint[[1]]; SE_joint <- joint[[2]]
  if(length(SE_joint)<length(est_joint)){SE_joint <- rep(NA,length(est_joint))} ## if no convergence
  MSEP_joint <- mean((joint[[3]]-gendat$true_randeff)^2)
  
  return(list(est_naive=est_naive,SE_naive=SE_naive,MSEP_naive=MSEP_naive, ## Naive
              est_jno0s=est_jno0s,SE_jno0s=SE_jno0s,MSEP_jno0s=MSEP_jno0s, ## Joint - No 0s
              est_joint=est_joint,SE_joint=SE_joint,MSEP_joint=MSEP_joint, ## Joint Model
              N=nrow(df), ## total No. of units
              N0s=sum(df2$Nk==0))) ## No. of empty clusters

}



#############################
##  Fit marginal models
fit_mods_marg <- function(gendat,size_plus=0,ZIP=TRUE){
  df <- gendat$data; df2 <- gendat$data_clstlvl

  ## GEE-EXCH
  gee <- geeglm(Yki ~ X1ki + X2ki,id=IDki,data=df,family=binomial,corstr="exchangeable")
  est_gee <- gee$coef; SE_gee <- summary(gee)$coef[,2]
  if(length(SE_gee)<length(est_gee)){SE_gee <- rep(NA,length(est_gee))} ## if no convergence

  ## IEE
  iee <- geeglm(Yki ~ X1ki + X2ki,id=IDki,data=df,family=binomial,corstr="independence")
  est_iee <- iee$coef; SE_iee <- summary(iee)$coef[,2]
  if(length(SE_iee)<length(est_iee)){SE_iee <- rep(NA,length(est_iee))} ## if no convergence

  ## WGEE
  wgee <- suppressWarnings(geeglm(Yki ~ X1ki + X2ki,id=IDki,data=df,family=binomial,corstr="independence",weights=ICWki)) ## without suppress warnings, get warning about non integer successes in binomial, due to weights
  est_wgee <- wgee$coef; SE_wgee <- summary(wgee)$coef[,2]
  if(length(SE_wgee)<length(est_wgee)){SE_wgee <- rep(NA,length(est_wgee))} ## if no convergence

  ## JMMICS
  jmm <- JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=df$IDki,Z0k=cbind(1,df2$X1k),minNk=size_plus,ZIP=ZIP,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad = 60)
  est_jmm <- jmm[[1]]; SE_jmm <- jmm[[2]]
  if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## if no convergence
  MSEP_jmm <- mean((jmm[[3]]-gendat$true_randeff)^2)
  
  return(list(est_gee=est_gee,SE_gee=SE_gee, ## exchangeable GEE
              est_iee=est_iee,SE_iee=SE_iee, ## IEE
              est_wgee=est_wgee,SE_wgee=SE_wgee, ## WGEE
              est_jmm=est_jmm,SE_jmm=SE_jmm,MSEP_jmm=MSEP_jmm, ## JMM
              N=nrow(df), ## total No. of units
              N0s=sum(df2$Nk==0))) ## No. of empty clusters
}


###########################
## Generate data and fit models
simICSfit <- function(gam=0, ## scaling factor in volume model 
                      m=0, ## minimum cluster size (0 or 1) ## MUST SET TO 0 FOR ZIP model
                      iter=0){ ## iteration number used to set seed for reproducibility
  
  if(iter!=0){ 
    set.seed(1234+10*iter) 
  } ## set seed
  
  ## generate conditional and marginal datasets
  cdat <- genICS(gamma_bk=gam,min_Nk=m,cond=TRUE,ZIP=TRUE,   K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma) ## conditional
  mdat <- genICS(gamma_bk=gam,min_Nk=m,cond=FALSE,ZIP=TRUE,  K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma) ## marginal
  
  ## fit conditional and marginal models
  cfit <- fit_mods_cond(cdat,m) 
  mfit <- fit_mods_marg(mdat,m)
  
  ## compile results
  res <- matrix(c(iter,
                  cfit$est_naive,cfit$SE_naive,cfit$MSEP_naive,
                  cfit$est_jno0s,cfit$SE_jno0s,cfit$MSEP_jno0s,
                  cfit$est_joint,cfit$SE_joint,cfit$MSEP_joint,
                  cfit$N,cfit$N0s,
                  mfit$est_gee,mfit$SE_gee,
                  mfit$est_iee,mfit$SE_iee,
                  mfit$est_wgee,mfit$SE_wgee,
                  mfit$est_jmm,mfit$SE_jmm,mfit$MSEP_jmm,
                  mfit$N,mfit$N0s),
                nrow=1)
  
  if(ncol(res)!=outputsize){ return(c(iter,rep(NA,outputsize-1)))} ## if incorrect number of results
  
  ## name columns
  colnames(res) <- c("*iter",
                     "B0est_naive","B1est_naive","B2est_naive","Sigest_naive","B0SE_naive","B1SE_naive","B2SE_naive","SigSE_naive","MSEP_naive",
                     "A0est_jno0s","A1est_jno0s","Gest_jno0s","B0est_jno0s","B1est_jno0s","B2est_jno0s","Sigest_jno0s","A0SE_jno0s","A1SE_jno0s","GSE_jno0s","B0SE_jno0s","B1SE_jno0s","B2SE_jno0s","SigSE_jno0s","MSEP_jno0s",
                     "A0est_joint","A1est_joint","Gest_joint","B0est_joint","B1est_joint","B2est_joint","Sigest_joint","E0est_joint","E1est_joint","A0SE_joint","A1SE_joint","GSE_joint","B0SE_joint","B1SE_joint","B2SE_joint","SigSE_joint","E0SE_joint","E1SE_joint","MSEP_joint",
                     "*N_cond","*N0s_cond",
                     "B0est_gee","B1est_gee","B2est_gee","B0SE_gee","B1SE_gee","B2SE_gee",
                     "B0est_iee","B1est_iee","B2est_iee","B0SE_iee","B1SE_iee","B2SE_iee",
                     "B0est_wgee","B1est_wgee","B2est_wgee","B0SE_wgee","B1SE_wgee","B2SE_wgee",
                     "A0est_jmm","A1est_jmm","Gest_jmm","B0est_jmm","B1est_jmm","B2est_jmm","Sigest_jmm","E0est_jmm","E1est_jmm","A0SE_jmm","A1SE_jmm","GSE_jmm","B0SE_jmm","B1SE_jmm","B2SE_jmm","SigSE_jmm","E0SE_jmm","E1SE_jmm","MSEP_jmm",
                     "*N_marg","*N0s_marg")
  
  ## reorder columns
  res <- res[,order(colnames(res))]
  
  return(res)
}











################
##    Loop    ##
################

## initialize results
# gam25 --> gamma= -0.25
# gamP25 --> gamma= +0.25
# m1 --> m= 1
#res_gam25_m1 <- res_gam10_m1 <- res_gam00_m1 <- c() ## res_gamP10_m1 <- res_gamP25_m1 <- c()
res_gam25_m0 <- res_gam10_m0 <- res_gam00_m0 <- c() ## res_gamP10_m0 <- res_gamP25_m0 <- c()


for (rr in loopvec){
  
  res_gam25_m0 <- rbind(res_gam25_m0,simICSfit(gam=-0.25,m=0,iter=rr))
  res_gam10_m0 <- rbind(res_gam10_m0,simICSfit(gam=-0.10,m=0,iter=rr))
  res_gam00_m0 <- rbind(res_gam00_m0,simICSfit(gam=0,m=0,iter=rr))
  # res_gamP10_m0 <- rbind(res_gamP10_m0,simICSfit(gam=0.10,m=0,iter=rr))
  # res_gamP25_m0 <- rbind(res_gamP25_m0,simICSfit(gam=0.25,m=0,iter=rr))
  
  print(rr)

}

write.table(res_gam25_m0,file=paste0(path,"results_ZIP_gam25_m0",suffix,".txt"))
write.table(res_gam10_m0,file=paste0(path,"results_ZIP_gam10_m0",suffix,".txt"))
write.table(res_gam00_m0,file=paste0(path,"results_ZIP_gam00_m0",suffix,".txt"))










# ################################
# ## gamma = -0.5    m = 1
# set.seed(1234+10*rr)
# cdat <- genICS(K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,gamma_bk=-0.5,sig.bk=sigma,min_Nk=1,cond=TRUE)
# mdat <- genICS(K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,gamma_bk=-0.5,sig.bk=sigma,min_Nk=1,cond=FALSE)
# cfit <- fit_mods_cond(cdat,1)
# mfit <- fit_mods_marg(mdat,1)

