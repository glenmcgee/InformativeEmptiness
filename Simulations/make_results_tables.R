###############################################################
##              Produce Results Tables and Plots             ##
###############################################################
## Edited: 11/21/2018
require(data.table)
require(xtable)
require(wesanderson)
require(tidyverse)
require(gridExtra)

## True parameter values
a0 <- 0.6; a1 <- -0.2 ## volume model coefficients
b0 <- -4; b1 <- 0.5; b2 <- 0.2 ## outcome model coefficients
e0 <- -2; e1 <- 1 ## zero-inflation model coefficients
#gamma <- 0,-0.1,-0.25 ## scaling factor for random intercept in volume model
## intercepts model
sigma <- 2 ## SD of rand int
## slopes model
sigma0 <- 2 ## SD of rand eff for X=0
sigma1 <- 1.5 ## SD of rand eff for X=1


## set path
path <- "./Simulations/Results/" ## path for results


## colors
wes_blue <- wes_palette(n=5, name="Darjeeling")[1]
wes_red <- wes_palette(n=5, name="Darjeeling")[3]
################################################
##                  Load Data                 ##
################################################
## Poisson
res_pois_gam25_m1 <- read.table(paste0(path,"results_poisson_gam25_m1.txt"),header = TRUE)
res_pois_gam10_m1 <- read.table(paste0(path,"results_poisson_gam10_m1.txt"),header = TRUE)
res_pois_gam00_m1 <- read.table(paste0(path,"results_poisson_gam00_m1.txt"),header = TRUE)
res_pois_gam25_m0 <- read.table(paste0(path,"results_poisson_gam25_m0.txt"),header = TRUE)
res_pois_gam10_m0 <- read.table(paste0(path,"results_poisson_gam10_m0.txt"),header = TRUE)
res_pois_gam00_m0 <- read.table(paste0(path,"results_poisson_gam00_m0.txt"),header = TRUE)
## Poisson + random slopes
res_pois_slopes_gam25_m1 <- read.table(paste0(path,"results_poisson_slopes_gam25_m1.txt"),header = TRUE)
res_pois_slopes_gam10_m1 <- read.table(paste0(path,"results_poisson_slopes_gam10_m1.txt"),header = TRUE)
res_pois_slopes_gam00_m1 <- read.table(paste0(path,"results_poisson_slopes_gam00_m1.txt"),header = TRUE)
res_pois_slopes_gam25_m0 <- read.table(paste0(path,"results_poisson_slopes_gam25_m0.txt"),header = TRUE)
res_pois_slopes_gam10_m0 <- read.table(paste0(path,"results_poisson_slopes_gam10_m0.txt"),header = TRUE)
res_pois_slopes_gam00_m0 <- read.table(paste0(path,"results_poisson_slopes_gam00_m0.txt"),header = TRUE)
## ZIP
res_ZIP_gam25_m0 <- read.table(paste0(path,"results_ZIP_gam25_m0.txt"),header = TRUE)
res_ZIP_gam10_m0 <- read.table(paste0(path,"results_ZIP_gam10_m0.txt"),header = TRUE)
res_ZIP_gam00_m0 <- read.table(paste0(path,"results_ZIP_gam00_m0.txt"),header = TRUE)
## ZIP +slopes
res_ZIP_slopes_gam25_m0 <- read.table(paste0(path,"results_ZIP_slopes_gam25_m0.txt"),header = TRUE)
res_ZIP_slopes_gam10_m0 <- read.table(paste0(path,"results_ZIP_slopes_gam10_m0.txt"),header = TRUE)
res_ZIP_slopes_gam00_m0 <- read.table(paste0(path,"results_ZIP_slopes_gam00_m0.txt"),header = TRUE)

################################################
##        Functions for Data Processing       ##
################################################
## function to identify outliers
get_outliers <- function(df,no_sd=4){
  ## get data
  df <- df[,order(colnames(df))] ## order data (just in case)
  ests <- df[,grep(pattern="est_",names(df))] ## get estimates
  meanmat <- matrix(apply(ests,2,mean,na.rm=T),byrow=T,nrow=nrow(ests),ncol=ncol(ests)) ## matrix of columnwise SDs
  sdmat <- matrix(apply(ests,2,sd,na.rm=T),byrow=T,nrow=nrow(ests),ncol=ncol(ests)) ## matrix of columnwise SDs
  # ## get parameters
  # true_vals <- get_truevals(gamma=gamma,ZIP=ZIP,slopes=slopes)
  # true_mat <- matrix(true_vals,nrow=nrow(ests),ncol=length(true_vals),byrow=TRUE) ## matrix of true parameters
  
  res <- abs(ests-meanmat)>no_sd*sdmat ## outliers
  which_outliers <- which(apply(res,1,sum)>0) ## which rows contain outliers
  return(which_outliers)
}

## function to get true parameters
get_truevals <- function(gamma=0,ZIP=FALSE,slopes=FALSE){
  
  ## add alphas and betas
  true_vals <- c(rep(c(a0,a1),each=3),  ## alpha only in joint models
                 rep(c(b0,b1,b2),each=7)) ## beta in all models
  
  ## add epsilons if ZIP model
  if(ZIP==TRUE){
    true_vals <- c(true_vals,
                   rep(c(e0,e1),each=2)) ## joint models EXCLUDING no0s (since it cant have zeros)
  }
  
  ## add gammas and sigmas as required
  if(slopes==FALSE){ ## random intercepts
    true_vals <- c(true_vals,
                   rep(gamma,each=3), ## only joint models
                   rep(sigma,each=4)) ## joint models AND naive (likelihood based methods)
    
  }else{ ## random slopes
    true_vals <- c(true_vals,
                   rep(c(gamma,gamma),each=3), ## only joint models ## currently gamma0=gamma1
                   rep(c(sigma0,sigma1),each=4)) ## joint models AND naive (likelihood based methods)
  }
  return(true_vals)
}




## function to compute pctbias and coverage (to later be averaged)
process_data <- function(df,gamma=0,ZIP=FALSE,slopes=FALSE){

  ## get data
  df <- df[,order(colnames(df))] ## order data (just in case)
  ests <- df[,grep(pattern="est_",names(df))] ## get estimates
  SEs <- df[,grep(pattern="SE_",names(df))] ## get standard errors
  
  ## get parameters
  true_vals <- get_truevals(gamma=gamma,ZIP=ZIP,slopes=slopes)
  true_mat <- matrix(true_vals,nrow=nrow(ests),ncol=length(true_vals),byrow=TRUE) ## matrix of true parameters
  
  ## get pctbias
  bias <- 100*(ests-true_mat)/true_mat   ## is true value covered by 95% Wald interval?
  colnames(bias) <- gsub(pattern="est",replacement="pctbias",colnames(bias))   ## rename columns
  
  ## get coveraage
  cvg <- 100*(true_mat>=ests-1.96*SEs & true_mat<=ests+1.96*SEs)   ## is true value covered by 95% Wald interval?
  colnames(cvg) <- gsub(pattern="est",replacement="cvg",colnames(cvg))   ## rename columns
  
  ## comine data
  df <- cbind(df,bias,cvg)
  
  return(df)
}


################################################
##              Data Processing               ##
################################################
## get and exclude outliers
# outliers_pois_gam25_m1 <- get_outliers(res_pois_gam25_m1); res_pois_gam25_m1 <- res_pois_gam25_m1[-outliers_pois_gam25_m1,]
# outliers_pois_gam10_m1 <- get_outliers(res_pois_gam10_m1); res_pois_gam10_m1 <- res_pois_gam10_m1[-outliers_pois_gam10_m1,]
# outliers_pois_gam00_m1 <- get_outliers(res_pois_gam00_m1); res_pois_gam00_m1 <- res_pois_gam00_m1[-outliers_pois_gam00_m1,]
# outliers_pois_gam25_m0 <- get_outliers(res_pois_gam25_m0); res_pois_gam25_m0 <- res_pois_gam25_m0[-outliers_pois_gam25_m0,]
# outliers_pois_gam10_m0 <- get_outliers(res_pois_gam10_m0); res_pois_gam10_m0 <- res_pois_gam10_m0[-outliers_pois_gam10_m0,]
# outliers_pois_gam00_m0 <- get_outliers(res_pois_gam00_m0); res_pois_gam00_m0 <- res_pois_gam00_m0[-outliers_pois_gam00_m0,]
# 
# outliers_pois_slopes_gam25_m1 <- get_outliers(res_pois_slopes_gam25_m1); res_pois_slopes_gam25_m1 <- res_pois_slopes_gam25_m1[-outliers_pois_slopes_gam25_m1,]
# outliers_pois_slopes_gam10_m1 <- get_outliers(res_pois_slopes_gam10_m1); res_pois_slopes_gam10_m1 <- res_pois_slopes_gam10_m1[-outliers_pois_slopes_gam10_m1,]
# outliers_pois_slopes_gam00_m1 <- get_outliers(res_pois_slopes_gam00_m1); res_pois_slopes_gam00_m1 <- res_pois_slopes_gam00_m1[-outliers_pois_slopes_gam00_m1,]
# outliers_pois_slopes_gam25_m0 <- get_outliers(res_pois_slopes_gam25_m0); res_pois_slopes_gam25_m0 <- res_pois_slopes_gam25_m0[-outliers_pois_slopes_gam25_m0,]
# outliers_pois_slopes_gam10_m0 <- get_outliers(res_pois_slopes_gam10_m0); res_pois_slopes_gam10_m0 <- res_pois_slopes_gam10_m0[-outliers_pois_slopes_gam10_m0,]
# outliers_pois_slopes_gam00_m0 <- get_outliers(res_pois_slopes_gam00_m0); res_pois_slopes_gam00_m0 <- res_pois_slopes_gam00_m0[-outliers_pois_slopes_gam00_m0,]
# 
# outliers_ZIP_gam25_m0 <- get_outliers(res_ZIP_gam25_m0); res_ZIP_gam25_m0 <- res_ZIP_gam25_m0[-outliers_ZIP_gam25_m0,]
# outliers_ZIP_gam10_m0 <- get_outliers(res_ZIP_gam10_m0); res_ZIP_gam10_m0 <- res_ZIP_gam10_m0[-outliers_ZIP_gam10_m0,]
# outliers_ZIP_gam00_m0 <- get_outliers(res_ZIP_gam00_m0); res_ZIP_gam00_m0 <- res_ZIP_gam00_m0[-outliers_ZIP_gam00_m0,]
# 
# outliers_ZIP_slopes_gam25_m0 <- get_outliers(res_ZIP_slopes_gam25_m0); res_ZIP_slopes_gam25_m0 <- res_ZIP_slopes_gam25_m0[-outliers_ZIP_slopes_gam25_m0,]
# outliers_ZIP_slopes_gam10_m0 <- get_outliers(res_ZIP_slopes_gam10_m0); res_ZIP_slopes_gam10_m0 <- res_ZIP_slopes_gam10_m0[-outliers_ZIP_slopes_gam10_m0,]
# outliers_ZIP_slopes_gam00_m0 <- get_outliers(res_ZIP_slopes_gam00_m0); res_ZIP_slopes_gam00_m0 <- res_ZIP_slopes_gam00_m0[-outliers_ZIP_slopes_gam00_m0,]



## Poisson
res_pois_gam00_m1 <- process_data(res_pois_gam00_m1,gamma=-0.00,ZIP=FALSE,slopes=FALSE)
res_pois_gam10_m1 <- process_data(res_pois_gam10_m1,gamma=-0.10,ZIP=FALSE,slopes=FALSE)
res_pois_gam25_m1 <- process_data(res_pois_gam25_m1,gamma=-0.25,ZIP=FALSE,slopes=FALSE)
res_pois_gam00_m0 <- process_data(res_pois_gam00_m0,gamma=-0.00,ZIP=FALSE,slopes=FALSE)
res_pois_gam10_m0 <- process_data(res_pois_gam10_m0,gamma=-0.10,ZIP=FALSE,slopes=FALSE)
res_pois_gam25_m0 <- process_data(res_pois_gam25_m0,gamma=-0.25,ZIP=FALSE,slopes=FALSE)
## Poisson+slopes
res_pois_slopes_gam00_m1 <- process_data(res_pois_slopes_gam00_m1,gamma=-0.00,ZIP=FALSE,slopes=TRUE)
res_pois_slopes_gam10_m1 <- process_data(res_pois_slopes_gam10_m1,gamma=-0.10,ZIP=FALSE,slopes=TRUE)
res_pois_slopes_gam25_m1 <- process_data(res_pois_slopes_gam25_m1,gamma=-0.25,ZIP=FALSE,slopes=TRUE)
res_pois_slopes_gam00_m0 <- process_data(res_pois_slopes_gam00_m0,gamma=-0.00,ZIP=FALSE,slopes=TRUE)
res_pois_slopes_gam10_m0 <- process_data(res_pois_slopes_gam10_m0,gamma=-0.10,ZIP=FALSE,slopes=TRUE)
res_pois_slopes_gam25_m0 <- process_data(res_pois_slopes_gam25_m0,gamma=-0.25,ZIP=FALSE,slopes=TRUE)
## ZIP
res_ZIP_gam00_m0 <- process_data(res_ZIP_gam00_m0,gamma=-0.00,ZIP=TRUE,slopes=FALSE)
res_ZIP_gam10_m0 <- process_data(res_ZIP_gam10_m0,gamma=-0.10,ZIP=TRUE,slopes=FALSE)
res_ZIP_gam25_m0 <- process_data(res_ZIP_gam25_m0,gamma=-0.25,ZIP=TRUE,slopes=FALSE)
## ZIP+slopes
res_ZIP_slopes_gam00_m0 <- process_data(res_ZIP_slopes_gam00_m0,gamma=-0.00,ZIP=TRUE,slopes=TRUE)
res_ZIP_slopes_gam10_m0 <- process_data(res_ZIP_slopes_gam10_m0,gamma=-0.10,ZIP=TRUE,slopes=TRUE)
res_ZIP_slopes_gam25_m0 <- process_data(res_ZIP_slopes_gam25_m0,gamma=-0.25,ZIP=TRUE,slopes=TRUE)


## list of data files
pois_ls <- list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1,
                res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0)
pois_slopes_ls <- list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1,
                       res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0)
ZIP_ls <- list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0)
ZIP_slopes_ls <- list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0)


################################################
##                   Tables                   ##
################################################


## handle NaNs
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}


make_tab <- function(df,stat="est",func=mean,ZIP=FALSE,slopes=FALSE){

  E0 <- E1 <- NULL
  A0 <- apply(cbind(NA,NA,NA,df[,paste0("A0",stat,"_jmm")],NA,df[,paste0("A0",stat,"_jno0s")],df[,paste0("A0",stat,"_joint")]),2,func,na.rm=TRUE)
  A1 <- apply(cbind(NA,NA,NA,df[,paste0("A1",stat,"_jmm")],NA,df[,paste0("A1",stat,"_jno0s")],df[,paste0("A1",stat,"_joint")]),2,func,na.rm=TRUE)
  B0 <- apply(cbind(df[,paste0("B0",stat,"_gee")],df[,paste0("B0",stat,"_iee")],df[,paste0("B0",stat,"_wgee")],df[,paste0("B0",stat,"_jmm")],df[,paste0("B0",stat,"_naive")],df[,paste0("B0",stat,"_jno0s")],df[,paste0("B0",stat,"_joint")]),2,func,na.rm=TRUE)
  B1 <- apply(cbind(df[,paste0("B1",stat,"_gee")],df[,paste0("B1",stat,"_iee")],df[,paste0("B1",stat,"_wgee")],df[,paste0("B1",stat,"_jmm")],df[,paste0("B1",stat,"_naive")],df[,paste0("B1",stat,"_jno0s")],df[,paste0("B1",stat,"_joint")]),2,func,na.rm=TRUE)
  B2 <- apply(cbind(df[,paste0("B2",stat,"_gee")],df[,paste0("B2",stat,"_iee")],df[,paste0("B2",stat,"_wgee")],df[,paste0("B2",stat,"_jmm")],df[,paste0("B2",stat,"_naive")],df[,paste0("B2",stat,"_jno0s")],df[,paste0("B2",stat,"_joint")]),2,func,na.rm=TRUE)
  Sig <- Sig0 <- Sig1 <- NULL
  G <- G0 <- G1 <- NULL
  
  if(slopes==FALSE){
    Sig <- apply(cbind(NA,NA,NA,df[,paste0("Sig",stat,"_jmm")],df[,paste0("Sig",stat,"_naive")],df[,paste0("Sig",stat,"_jno0s")],df[,paste0("Sig",stat,"_joint")]),2,func,na.rm=TRUE)
    G <- apply(cbind(NA,NA,NA,df[,paste0("G",stat,"_jmm")],NA,df[,paste0("G",stat,"_jno0s")],df[,paste0("G",stat,"_joint")]),2,func,na.rm=TRUE)
  }else{
    Sig0 <- apply(cbind(NA,NA,NA,df[,paste0("Sig0",stat,"_jmm")],df[,paste0("Sig0",stat,"_naive")],df[,paste0("Sig0",stat,"_jno0s")],df[,paste0("Sig0",stat,"_joint")]),2,func,na.rm=TRUE)
    Sig1 <- apply(cbind(NA,NA,NA,df[,paste0("Sig1",stat,"_jmm")],df[,paste0("Sig1",stat,"_naive")],df[,paste0("Sig1",stat,"_jno0s")],df[,paste0("Sig1",stat,"_joint")]),2,func,na.rm=TRUE)
    G0 <- apply(cbind(NA,NA,NA,df[,paste0("G0",stat,"_jmm")],NA,df[,paste0("G0",stat,"_jno0s")],df[,paste0("G0",stat,"_joint")]),2,func,na.rm=TRUE)
    G1 <- apply(cbind(NA,NA,NA,df[,paste0("G1",stat,"_jmm")],NA,df[,paste0("G1",stat,"_jno0s")],df[,paste0("G1",stat,"_joint")]),2,func,na.rm=TRUE)
  }

  if(ZIP==TRUE){
    E0 <-   apply(cbind(NA,NA,NA,df[,paste0("E0",stat,"_jmm")],NA,NA,df[,paste0("E0",stat,"_joint")]),2,func,na.rm=TRUE)
    E1 <-   apply(cbind(NA,NA,NA,df[,paste0("E1",stat,"_jmm")],NA,NA,df[,paste0("E1",stat,"_joint")]),2,func,na.rm=TRUE)
  }

  sum_df <- rbind(E0,E1,A0,A1,B0,B1,B2,Sig,Sig0,Sig1,G,G0,G1)
  colnames(sum_df) <- c("GEE","IEE","WGEE","JMM","Naive","Joint No 0s","Joint")
  sum_df[is.nan(sum_df)] <- NA
  
  return(data.frame(sum_df))
}

combine_tabs <- function(ls,gams,dig=2){

  for(ll in 1:length(ls)){
    ls[[ll]] <- round(ls[[ll]],dig) ## round
    ls[[ll]] <- data.frame(cbind(index=seq(1,nrow(ls[[ll]])),Param=rownames(ls[[ll]]),Gamma=gams[ll],ls[[ll]])) ## add index for reordering
  }
  res <- data.frame(rbindlist(ls)[order(index),-1]) ## combine
  rownames(res) <- paste(res$Param,res$Gamma) ## unique row name
  res <- data.frame(res) ## remove redundant column
  
  return(res)
}

summary_tab <- function(ls,gams=c(0,-0.1,-0.25),stat="est",func=mean,ZIP=FALSE,slopes=FALSE,dig=2){
  
  ## 
  for(ll in 1:length(ls)){
    ls[[ll]] <- make_tab(ls[[ll]],stat,func,ZIP,slopes) ## make individual tables
  }
  
  res <- combine_tabs(ls,gams,dig)
  
  return(res)
  
}

MSEP_tab <- function(ls,gams=c(0,-0.1,-0.25),dig=2){
  res <- c()
  for(ll in 1:length(ls)){
    df <- ls[[ll]]
    res <- rbind(res,round(c(median(df$MSEP_jmm),median(df$MSEP_naive),median(df$MSEP_jno0s),median(df$MSEP_joint)),dig))
  }
  colnames(res) <- c("JMM","Naive","Joint No 0s","Joint")
  res <- cbind(gams,res)
  return(res)
}


## mean estimates
mean_pois_m1 <- summary_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1),c(0,-0.1,-0.25),"est",mean,ZIP=FALSE,slopes=FALSE,dig=2)
mean_pois_slopes_m1 <- summary_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1),c(0,-0.1,-0.25),"est",mean,ZIP=FALSE,slopes=TRUE,dig=2)
mean_pois_m0 <- summary_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0),c(0,-0.1,-0.25),"est",mean,ZIP=FALSE,slopes=FALSE,dig=2)
mean_pois_slopes_m0 <- summary_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0),c(0,-0.1,-0.25),"est",mean,ZIP=FALSE,slopes=TRUE,dig=2)
mean_ZIP_m0 <- summary_tab(list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0),c(0,-0.1,-0.25),"est",mean,ZIP=TRUE,slopes=FALSE,dig=2)
mean_ZIP_slopes_m0 <- summary_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0),c(0,-0.1,-0.25),"est",mean,ZIP=TRUE,slopes=TRUE,dig=2)

## median estimates
median_pois_m1 <- summary_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1),c(0,-0.1,-0.25),"est",median,ZIP=FALSE,slopes=FALSE,dig=2)
median_pois_slopes_m1 <- summary_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1),c(0,-0.1,-0.25),"est",median,ZIP=FALSE,slopes=TRUE,dig=2)
median_pois_m0 <- summary_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0),c(0,-0.1,-0.25),"est",median,ZIP=FALSE,slopes=FALSE,dig=2)
median_pois_slopes_m0 <- summary_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0),c(0,-0.1,-0.25),"est",median,ZIP=FALSE,slopes=TRUE,dig=2)
median_ZIP_m0 <- summary_tab(list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0),c(0,-0.1,-0.25),"est",median,ZIP=TRUE,slopes=FALSE,dig=2)
median_ZIP_slopes_m0 <- summary_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0),c(0,-0.1,-0.25),"est",median,ZIP=TRUE,slopes=TRUE,dig=2)

## pct bias
pctbias_pois_m1 <- summary_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1),c(0,-0.1,-0.25),"pctbias",mean,ZIP=FALSE,slopes=FALSE,dig=0)
pctbias_pois_slopes_m1 <- summary_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1),c(0,-0.1,-0.25),"pctbias",mean,ZIP=FALSE,slopes=TRUE,dig=0)
pctbias_pois_m0 <- summary_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0),c(0,-0.1,-0.25),"pctbias",mean,ZIP=FALSE,slopes=FALSE,dig=0)
pctbias_pois_slopes_m0 <- summary_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0),c(0,-0.1,-0.25),"pctbias",mean,ZIP=FALSE,slopes=TRUE,dig=0)
pctbias_ZIP_m0 <- summary_tab(list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0),c(0,-0.1,-0.25),"pctbias",mean,ZIP=TRUE,slopes=FALSE,dig=0)
pctbias_ZIP_slopes_m0 <- summary_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0),c(0,-0.1,-0.25),"pctbias",mean,ZIP=TRUE,slopes=TRUE,dig=0)

## median pct bias
medpctbias_pois_m1 <- summary_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1),c(0,-0.1,-0.25),"pctbias",median,ZIP=FALSE,slopes=FALSE,dig=0)
medpctbias_pois_slopes_m1 <- summary_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1),c(0,-0.1,-0.25),"pctbias",median,ZIP=FALSE,slopes=TRUE,dig=0)
medpctbias_pois_m0 <- summary_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0),c(0,-0.1,-0.25),"pctbias",median,ZIP=FALSE,slopes=FALSE,dig=0)
medpctbias_pois_slopes_m0 <- summary_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0),c(0,-0.1,-0.25),"pctbias",median,ZIP=FALSE,slopes=TRUE,dig=0)
medpctbias_ZIP_m0 <- summary_tab(list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0),c(0,-0.1,-0.25),"pctbias",median,ZIP=TRUE,slopes=FALSE,dig=0)
medpctbias_ZIP_slopes_m0 <- summary_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0),c(0,-0.1,-0.25),"pctbias",median,ZIP=TRUE,slopes=TRUE,dig=0)

## 95% coverage
cvg_pois_m1 <- summary_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1),c(0,-0.1,-0.25),"cvg",mean,ZIP=FALSE,slopes=FALSE,dig=0)
cvg_pois_slopes_m1 <- summary_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1),c(0,-0.1,-0.25),"cvg",mean,ZIP=FALSE,slopes=TRUE,dig=0)
cvg_pois_m0 <- summary_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0),c(0,-0.1,-0.25),"cvg",mean,ZIP=FALSE,slopes=FALSE,dig=0)
cvg_pois_slopes_m0 <- summary_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0),c(0,-0.1,-0.25),"cvg",mean,ZIP=FALSE,slopes=TRUE,dig=0)
cvg_ZIP_m0 <- summary_tab(list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0),c(0,-0.1,-0.25),"cvg",mean,ZIP=TRUE,slopes=FALSE,dig=0)
cvg_ZIP_slopes_m0 <- summary_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0),c(0,-0.1,-0.25),"cvg",mean,ZIP=TRUE,slopes=TRUE,dig=0)

## SD of estimates
sd_pois_m1 <- summary_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1),c(0,-0.1,-0.25),"est",sd,ZIP=FALSE,slopes=FALSE,dig=2)
sd_pois_slopes_m1 <- summary_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1),c(0,-0.1,-0.25),"est",sd,ZIP=FALSE,slopes=TRUE,dig=2)
sd_pois_m0 <- summary_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0),c(0,-0.1,-0.25),"est",sd,ZIP=FALSE,slopes=FALSE,dig=2)
sd_pois_slopes_m0 <- summary_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0),c(0,-0.1,-0.25),"est",sd,ZIP=FALSE,slopes=TRUE,dig=2)
sd_ZIP_m0 <- summary_tab(list(res_ZIP_gam00_m0,res_ZIP_gam10_m0,res_ZIP_gam25_m0),c(0,-0.1,-0.25),"est",sd,ZIP=TRUE,slopes=FALSE,dig=2)
sd_ZIP_slopes_m0 <- summary_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0),c(0,-0.1,-0.25),"est",sd,ZIP=TRUE,slopes=TRUE,dig=2)



########################################
## JUST BETA TABLE (PCT BIAS)
marg_id <- 3:6 ## which cols for marginal
cond_id <- 7:9

## marginal bias table
Beta_tab_marg <- rbind(cbind(medpctbias_pois_m1[grep(pattern="B",rownames(medpctbias_pois_m1)),c(1:2,marg_id)],medpctbias_pois_m0[grep(pattern="B",rownames(medpctbias_pois_m0)),c(marg_id)]),
                       cbind(medpctbias_pois_slopes_m1[grep(pattern="B",rownames(medpctbias_pois_slopes_m1)),c(1:2,marg_id)],medpctbias_pois_slopes_m0[grep(pattern="B",rownames(medpctbias_pois_slopes_m0)),c(marg_id)]) )
Beta_tab_marg_ZIP <- rbind(cbind(medpctbias_ZIP_m0[grep(pattern="B",rownames(medpctbias_ZIP_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(marg_id)),medpctbias_ZIP_m0[grep(pattern="B",rownames(medpctbias_ZIP_m0)),c(marg_id)]),
                           cbind(medpctbias_ZIP_slopes_m0[grep(pattern="B",rownames(medpctbias_ZIP_slopes_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(marg_id)),medpctbias_ZIP_slopes_m0[grep(pattern="B",rownames(medpctbias_ZIP_slopes_m0)),c(marg_id)]))
colnames(Beta_tab_marg_ZIP) <- colnames(Beta_tab_marg)
Beta_tab_marg <- rbind(Beta_tab_marg,Beta_tab_marg_ZIP)

## conditional bias table
Beta_tab_cond <- rbind(cbind(medpctbias_pois_m1[grep(pattern="B",rownames(medpctbias_pois_m1)),c(1:2,cond_id)],medpctbias_pois_m0[grep(pattern="B",rownames(medpctbias_pois_m0)),c(cond_id)]),
                       cbind(medpctbias_pois_slopes_m1[grep(pattern="B",rownames(medpctbias_pois_slopes_m1)),c(1:2,cond_id)],medpctbias_pois_slopes_m0[grep(pattern="B",rownames(medpctbias_pois_slopes_m0)),c(cond_id)]) )
Beta_tab_cond_ZIP <- rbind(cbind(medpctbias_ZIP_m0[grep(pattern="B",rownames(medpctbias_ZIP_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(cond_id)),medpctbias_ZIP_m0[grep(pattern="B",rownames(medpctbias_ZIP_m0)),c(cond_id)]),
                           cbind(medpctbias_ZIP_slopes_m0[grep(pattern="B",rownames(medpctbias_ZIP_slopes_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(cond_id)),medpctbias_ZIP_slopes_m0[grep(pattern="B",rownames(medpctbias_ZIP_slopes_m0)),c(cond_id)]))
colnames(Beta_tab_cond_ZIP) <- colnames(Beta_tab_cond)
Beta_tab_cond <- rbind(Beta_tab_cond,Beta_tab_cond_ZIP)

## marginal table
medBeta_tab_marg <- rbind(cbind(median_pois_m1[grep(pattern="B",rownames(median_pois_m1)),c(1:2,marg_id)],median_pois_m0[grep(pattern="B",rownames(median_pois_m0)),c(marg_id)]),
                       cbind(median_pois_slopes_m1[grep(pattern="B",rownames(median_pois_slopes_m1)),c(1:2,marg_id)],median_pois_slopes_m0[grep(pattern="B",rownames(median_pois_slopes_m0)),c(marg_id)]) )
medBeta_tab_marg_ZIP <- rbind(cbind(median_ZIP_m0[grep(pattern="B",rownames(median_ZIP_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(marg_id)),median_ZIP_m0[grep(pattern="B",rownames(median_ZIP_m0)),c(marg_id)]),
                           cbind(median_ZIP_slopes_m0[grep(pattern="B",rownames(median_ZIP_slopes_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(marg_id)),median_ZIP_slopes_m0[grep(pattern="B",rownames(median_ZIP_slopes_m0)),c(marg_id)]))
colnames(medBeta_tab_marg_ZIP) <- colnames(medBeta_tab_marg)
medBeta_tab_marg <- rbind(medBeta_tab_marg,medBeta_tab_marg_ZIP)

## conditional table
medBeta_tab_cond <- rbind(cbind(median_pois_m1[grep(pattern="B",rownames(median_pois_m1)),c(1:2,cond_id)],median_pois_m0[grep(pattern="B",rownames(median_pois_m0)),c(cond_id)]),
                       cbind(median_pois_slopes_m1[grep(pattern="B",rownames(median_pois_slopes_m1)),c(1:2,cond_id)],median_pois_slopes_m0[grep(pattern="B",rownames(median_pois_slopes_m0)),c(cond_id)]) )
medBeta_tab_cond_ZIP <- rbind(cbind(median_ZIP_m0[grep(pattern="B",rownames(median_ZIP_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(cond_id)),median_ZIP_m0[grep(pattern="B",rownames(median_ZIP_m0)),c(cond_id)]),
                           cbind(median_ZIP_slopes_m0[grep(pattern="B",rownames(median_ZIP_slopes_m0)),c(1:2)],matrix(NA,nrow=9,ncol=length(cond_id)),median_ZIP_slopes_m0[grep(pattern="B",rownames(median_ZIP_slopes_m0)),c(cond_id)]))
colnames(medBeta_tab_cond_ZIP) <- colnames(medBeta_tab_cond)
medBeta_tab_cond <- rbind(medBeta_tab_cond,medBeta_tab_cond_ZIP)


########################################
## MSEP 
MSEP_pois_m1 <- MSEP_tab(list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1))
MSEP_pois_m0 <- MSEP_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0))
MSEP_pois_slopes_m1 <- MSEP_tab(list(res_pois_slopes_gam00_m1,res_pois_slopes_gam10_m1,res_pois_slopes_gam25_m1))
MSEP_pois_slopes_m0 <- MSEP_tab(list(res_pois_slopes_gam00_m0,res_pois_slopes_gam10_m0,res_pois_slopes_gam25_m0))
MSEP_ZIP <- MSEP_tab(list(res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0))
MSEP_ZIP_slopes <- MSEP_tab(list(res_ZIP_slopes_gam00_m0,res_ZIP_slopes_gam10_m0,res_ZIP_slopes_gam25_m0))

MSEP_table <- rbind(cbind(MSEP_pois_m1,MSEP_pois_m0[,-1]),
                    cbind(MSEP_pois_slopes_m1,MSEP_pois_slopes_m0[,-1]),
                    cbind(MSEP_ZIP[,1],matrix(NA,nrow=3,ncol=4),MSEP_ZIP[,-1]),
                    cbind(MSEP_ZIP[,1],matrix(NA,nrow=3,ncol=4),MSEP_ZIP_slopes[,-1]))
MSEP_table <- data.frame(cbind(rep(c("Poisson","Poisson Slopes","ZIP","ZIP Slopes"),each=3),MSEP_table))




################################################
##                  Boxplots                  ##
################################################
make_box <- function(df,caption=""){#,stat="est",func=mean,ZIP=FALSE,slopes=FALSE){
  
  B0 <- df[,grep(pattern="B0est_",names(df))]
  colnames(B0) <- c("1GEE","2IEE","4JMM","6No 0s","7Joint","5Naive","3WGEE")
  B0 <- B0%>%gather("Group","Est",1:7)
  B0$color <- wes_red; B0$color[B0$Group %in% c("1GEE","2IEE","3WGEE","4JMM")] <- wes_blue
  B0$opacity <- 0.75; B0$opacity[B0$Group %in% c("4JMM","7Joint")] <- 1
  
  B1 <- df[,grep(pattern="B1est_",names(df))]
  colnames(B1) <- c("1GEE","2IEE","4JMM","6No 0s","7Joint","5Naive","3WGEE")
  B1 <- B1%>%gather("Group","Est",1:7)
  B1$color <- wes_red; B1$color[B1$Group %in% c("1GEE","2IEE","3WGEE","4JMM")] <- wes_blue
  B1$opacity <- 0.75; B1$opacity[B1$Group %in% c("4JMM","7Joint")] <- 1
  
  B2 <- df[,grep(pattern="B2est_",names(df))]
  colnames(B2) <- c("1GEE","2IEE","4JMM","6No 0s","7Joint","5Naive","3WGEE")
  B2 <- B2%>%gather("Group","Est",1:7)
  B2$color <- wes_red; B2$color[B2$Group %in% c("1GEE","2IEE","3WGEE","4JMM")] <- wes_blue
  B2$opacity <- 0.75; B2$opacity[B2$Group %in% c("4JMM","7Joint")] <- 1
  
  
  
  box_B0 <- ggplot(B0,aes(factor(Group),Est,fill=color,alpha=as.factor(opacity)))+
    geom_boxplot(aes())+
    geom_hline(yintercept=b0)+
    scale_y_continuous(limits=c(b0-2,b0+2))+
    scale_fill_discrete(guide=F)+
    scale_alpha_discrete(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    theme_classic() +
    ggtitle(paste(caption)) +
    theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 10, b = -20)))

  
  box_B1 <- ggplot(B1,aes(factor(Group),Est,fill=color,alpha=as.factor(opacity)))+
    geom_boxplot(aes())+
    geom_hline(yintercept=b1)+
    scale_y_continuous(limits=c(b1-1.5,b1+1.5))+
    scale_fill_discrete(guide=F)+
    scale_alpha_discrete(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    theme_classic()+
    ggtitle(paste(caption)) +
    theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 10, b = -20)))

  box_B2 <- ggplot(B2,aes(factor(Group),Est,fill=color,alpha=as.factor(opacity)))+
    geom_boxplot(aes())+
    geom_hline(yintercept=b2)+
    scale_y_continuous(limits=c(b2-1,b2+1))+
    scale_fill_discrete(guide=F)+
    scale_alpha_discrete(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    theme_classic()+
    ggtitle(paste(caption)) +
    theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 10, b = -20)))
  
  return(list(B0=box_B0,B1=box_B1,B2=box_B2))
}


box_pois_gam00_m1 <- make_box(res_pois_gam00_m1,"gamma=0, m=1")
box_pois_gam10_m1 <- make_box(res_pois_gam10_m1,"gamma=-0.1, m=1")
box_pois_gam25_m1 <- make_box(res_pois_gam25_m1,"gamma=-0.25, m=1")
box_pois_gam00_m0 <- make_box(res_pois_gam00_m0,"gamma=0, m=0")
box_pois_gam10_m0 <- make_box(res_pois_gam10_m0,"gamma=-0.1, m=0")
box_pois_gam25_m0 <- make_box(res_pois_gam25_m0,"gamma=-0.25, m=0")

box_pois_slopes_gam00_m1 <- make_box(res_pois_slopes_gam00_m1,"gamma=0, m=1")
box_pois_slopes_gam10_m1 <- make_box(res_pois_slopes_gam10_m1,"gamma=-0.1, m=1")
box_pois_slopes_gam25_m1 <- make_box(res_pois_slopes_gam25_m1,"gamma=-0.25, m=1")
box_pois_slopes_gam00_m0 <- make_box(res_pois_slopes_gam00_m0,"gamma=0, m=0")
box_pois_slopes_gam10_m0 <- make_box(res_pois_slopes_gam10_m0,"gamma=-0.1, m=0")
box_pois_slopes_gam25_m0 <- make_box(res_pois_slopes_gam25_m0,"gamma=-0.25, m=0")


box_ZIP_gam00_m0 <- make_box(res_ZIP_gam00_m0,"gamma=0, ZIP")
box_ZIP_gam10_m0 <- make_box(res_ZIP_gam10_m0,"gamma=-0.1, ZIP")
box_ZIP_gam25_m0 <- make_box(res_ZIP_gam25_m0,"gamma=-0.25, ZIP")

box_ZIP_slopes_gam00_m0 <- make_box(res_ZIP_slopes_gam00_m0,"gamma=0, ZIP")
box_ZIP_slopes_gam10_m0 <- make_box(res_ZIP_slopes_gam10_m0,"gamma=-0.1, ZIP")
box_ZIP_slopes_gam25_m0 <- make_box(res_ZIP_slopes_gam25_m0,"gamma=-0.25, ZIP")



################################################
##                 Nk and 0s                  ##
################################################
## cond and marg
meanNk_pois_gam00_m1 <- c(mean(res_pois_gam00_m1$X.N_cond),mean(res_pois_gam00_m1$X.N_marg))
meanNk_pois_gam10_m1 <- c(mean(res_pois_gam10_m1$X.N_cond),mean(res_pois_gam10_m1$X.N_marg))
meanNk_pois_gam25_m1 <- c(mean(res_pois_gam25_m1$X.N_cond),mean(res_pois_gam25_m1$X.N_marg))
meanNk_pois_gam00_m0 <- c(mean(res_pois_gam00_m0$X.N_cond),mean(res_pois_gam00_m0$X.N_marg))
meanNk_pois_gam10_m0 <- c(mean(res_pois_gam10_m0$X.N_cond),mean(res_pois_gam10_m0$X.N_marg))
meanNk_pois_gam25_m0 <- c(mean(res_pois_gam25_m0$X.N_cond),mean(res_pois_gam25_m0$X.N_marg))

meanNk_pois_slopes_gam00_m1 <- c(mean(res_pois_slopes_gam00_m1$X.N_cond),mean(res_pois_slopes_gam00_m1$X.N_marg))
meanNk_pois_slopes_gam10_m1 <- c(mean(res_pois_slopes_gam10_m1$X.N_cond),mean(res_pois_slopes_gam10_m1$X.N_marg))
meanNk_pois_slopes_gam25_m1 <- c(mean(res_pois_slopes_gam25_m1$X.N_cond),mean(res_pois_slopes_gam25_m1$X.N_marg))
meanNk_pois_slopes_gam00_m0 <- c(mean(res_pois_slopes_gam00_m0$X.N_cond),mean(res_pois_slopes_gam00_m0$X.N_marg))
meanNk_pois_slopes_gam10_m0 <- c(mean(res_pois_slopes_gam10_m0$X.N_cond),mean(res_pois_slopes_gam10_m0$X.N_marg))
meanNk_pois_slopes_gam25_m0 <- c(mean(res_pois_slopes_gam25_m0$X.N_cond),mean(res_pois_slopes_gam25_m0$X.N_marg))

meanNk_ZIP_gam00_m0 <- c(mean(res_ZIP_gam00_m0$X.N_cond),mean(res_ZIP_gam00_m0$X.N_marg))
meanNk_ZIP_gam10_m0 <- c(mean(res_ZIP_gam10_m0$X.N_cond),mean(res_ZIP_gam10_m0$X.N_marg))
meanNk_ZIP_gam25_m0 <- c(mean(res_ZIP_gam25_m0$X.N_cond),mean(res_ZIP_gam25_m0$X.N_marg))

meanNk_ZIP_slopes_gam00_m0 <- c(mean(res_ZIP_slopes_gam00_m0$X.N_cond),mean(res_ZIP_slopes_gam00_m0$X.N_marg))
meanNk_ZIP_slopes_gam10_m0 <- c(mean(res_ZIP_slopes_gam10_m0$X.N_cond),mean(res_ZIP_slopes_gam10_m0$X.N_marg))
meanNk_ZIP_slopes_gam25_m0 <- c(mean(res_ZIP_slopes_gam25_m0$X.N_cond),mean(res_ZIP_slopes_gam25_m0$X.N_marg))

## make table of sizes
meanNks <- rbind(cbind(rbind(meanNk_pois_gam00_m1,meanNk_pois_gam10_m1,meanNk_pois_gam25_m1),
                         rbind(meanNk_pois_gam00_m0,meanNk_pois_gam10_m0,meanNk_pois_gam25_m0)),
                   
                   cbind(rbind(meanNk_pois_slopes_gam00_m1,meanNk_pois_slopes_gam10_m1,meanNk_pois_slopes_gam25_m1),
                         rbind(meanNk_pois_slopes_gam00_m0,meanNk_pois_slopes_gam10_m0,meanNk_pois_slopes_gam25_m0)),
                   
                 cbind(matrix(NA,nrow=3,ncol=2),rbind(meanNk_ZIP_gam00_m0,meanNk_ZIP_gam10_m0,meanNk_ZIP_gam25_m0)),
                 cbind(matrix(NA,nrow=3,ncol=2),rbind(meanNk_ZIP_slopes_gam00_m0,meanNk_ZIP_slopes_gam10_m0,meanNk_ZIP_slopes_gam25_m0)))
meanNks <- round(meanNks)
colnames(meanNks) <- c("Cond","Marg","Cond0","Marg1")
meanNks <- cbind(models=rep(c("Poisson","Poisson Slopes","ZIP","ZIP Slopes"),each=3),gamma=rep(c(0,-0.10,-0.25),4),meanNks)

## cond and marg
pctzero_pois_gam00_m0 <- 100*(1/nrow(res_pois_gam00_m0))*c(mean(res_pois_gam00_m0$X.N0s_cond),mean(res_pois_gam00_m0$X.N0s_marg))
pctzero_pois_gam10_m0 <- 100*(1/nrow(res_pois_gam10_m0))*c(mean(res_pois_gam10_m0$X.N0s_cond),mean(res_pois_gam10_m0$X.N0s_marg))
pctzero_pois_gam25_m0 <- 100*(1/nrow(res_pois_gam25_m0))*c(mean(res_pois_gam25_m0$X.N0s_cond),mean(res_pois_gam25_m0$X.N0s_marg))

pctzero_pois_slopes_gam00_m0 <- 100*(1/nrow(res_pois_slopes_gam00_m0))*c(mean(res_pois_slopes_gam00_m0$X.N0s_cond),mean(res_pois_slopes_gam00_m0$X.N0s_marg))
pctzero_pois_slopes_gam10_m0 <- 100*(1/nrow(res_pois_slopes_gam10_m0))*c(mean(res_pois_slopes_gam10_m0$X.N0s_cond),mean(res_pois_slopes_gam10_m0$X.N0s_marg))
pctzero_pois_slopes_gam25_m0 <- 100*(1/nrow(res_pois_slopes_gam25_m0))*c(mean(res_pois_slopes_gam25_m0$X.N0s_cond),mean(res_pois_slopes_gam25_m0$X.N0s_marg))

pctzero_ZIP_gam00_m0 <- 100*(1/nrow(res_ZIP_gam00_m0))*c(mean(res_ZIP_gam00_m0$X.N0s_cond),mean(res_ZIP_gam00_m0$X.N0s_marg))
pctzero_ZIP_gam10_m0 <- 100*(1/nrow(res_ZIP_gam10_m0))*c(mean(res_ZIP_gam10_m0$X.N0s_cond),mean(res_ZIP_gam10_m0$X.N0s_marg))
pctzero_ZIP_gam25_m0 <- 100*(1/nrow(res_ZIP_gam25_m0))*c(mean(res_ZIP_gam25_m0$X.N0s_cond),mean(res_ZIP_gam25_m0$X.N0s_marg))

pctzero_ZIP_slopes_gam00_m0 <- 100*(1/nrow(res_ZIP_slopes_gam00_m0))*c(mean(res_ZIP_slopes_gam00_m0$X.N0s_cond),mean(res_ZIP_slopes_gam00_m0$X.N0s_marg))
pctzero_ZIP_slopes_gam10_m0 <- 100*(1/nrow(res_ZIP_slopes_gam10_m0))*c(mean(res_ZIP_slopes_gam10_m0$X.N0s_cond),mean(res_ZIP_slopes_gam10_m0$X.N0s_marg))
pctzero_ZIP_slopes_gam25_m0 <- 100*(1/nrow(res_ZIP_slopes_gam25_m0))*c(mean(res_ZIP_slopes_gam25_m0$X.N0s_cond),mean(res_ZIP_slopes_gam25_m0$X.N0s_marg))

pctzeros_tab <- rbind(pctzero_pois_gam00_m0,pctzero_pois_gam10_m0,pctzero_pois_gam25_m0,
                      pctzero_pois_slopes_gam00_m0,pctzero_pois_slopes_gam10_m0,pctzero_pois_slopes_gam25_m0,
                      pctzero_ZIP_gam00_m0,pctzero_ZIP_gam10_m0,pctzero_ZIP_gam25_m0,
                      pctzero_ZIP_slopes_gam00_m0,pctzero_ZIP_slopes_gam10_m0,pctzero_ZIP_slopes_gam25_m0)
pctzeros_tab <- round(pctzeros_tab)
colnames(pctzeros_tab) <- c("Cond","Marg")
pctzeros_tab <- cbind(models=rep(c("Poisson","Poisson Slopes","ZIP","ZIP Slopes"),each=3),gamma=rep(c(0,-0.10,-0.25),4),pctzeros_tab)


## function for size hists
make_sizehist <- function(df,caption=""){

  marg <- ggplot(df,aes(X.N_marg))+
    geom_histogram(bins=60)+
    scale_x_continuous(limits=c(3000,6000))+
    ylab("No of populations")+xlab("Population size")+
    theme_classic() +
    ggtitle(paste("Marginal -",caption))
  
  cond <- ggplot(df,aes(X.N_cond))+
    geom_histogram(bins=60)+
    scale_x_continuous(limits=c(3000,6000))+
    ylab("No of populations")+xlab("Population size")+
    theme_classic() +
    ggtitle(paste("Conditional -",caption))
  
  return(list(marg=marg,cond=cond))
  
}

size_hist_pois_gam00_m1 <- make_sizehist(res_pois_gam00_m1,"gamma=0 - m=1")
size_hist_pois_gam10_m1 <- make_sizehist(res_pois_gam10_m1,"gamma=-0.1 - m=1")
size_hist_pois_gam25_m1 <- make_sizehist(res_pois_gam25_m1,"gamma=-0.25 - m=1")
size_hist_pois_gam00_m0 <- make_sizehist(res_pois_gam00_m0,"gamma=0 - m=0")
size_hist_pois_gam10_m0 <- make_sizehist(res_pois_gam10_m0,"gamma=-0.1 - m=0")
size_hist_pois_gam25_m0 <- make_sizehist(res_pois_gam25_m0,"gamma=-0.25 - m=0")

size_hist_pois_slopes_gam00_m1 <- make_sizehist(res_pois_slopes_gam00_m1,"gamma=0 - m=1")
size_hist_pois_slopes_gam10_m1 <- make_sizehist(res_pois_slopes_gam10_m1,"gamma=-0.1 - m=1")
size_hist_pois_slopes_gam25_m1 <- make_sizehist(res_pois_slopes_gam25_m1,"gamma=-0.25 - m=1")
size_hist_pois_slopes_gam00_m0 <- make_sizehist(res_pois_slopes_gam00_m0,"gamma=0 - m=0")
size_hist_pois_slopes_gam10_m0 <- make_sizehist(res_pois_slopes_gam10_m0,"gamma=-0.1 - m=0")
size_hist_pois_slopes_gam25_m0 <- make_sizehist(res_pois_slopes_gam25_m0,"gamma=-0.25 - m=0")

size_hist_ZIP_gam00_m0 <- make_sizehist(res_ZIP_gam00_m0,"gamma=0 - m=0")
size_hist_ZIP_gam10_m0 <- make_sizehist(res_ZIP_gam10_m0,"gamma=-0.1 - m=0")
size_hist_ZIP_gam25_m0 <- make_sizehist(res_ZIP_gam25_m0,"gamma=-0.25 - m=0")

size_hist_ZIP_slopes_gam00_m0 <- make_sizehist(res_ZIP_slopes_gam00_m0,"gamma=0 - m=0")
size_hist_ZIP_slopes_gam10_m0 <- make_sizehist(res_ZIP_slopes_gam10_m0,"gamma=-0.1 - m=0")
size_hist_ZIP_slopes_gam25_m0 <- make_sizehist(res_ZIP_slopes_gam25_m0,"gamma=-0.25 - m=0")


## function for size hists
make_zerohist <- function(df,caption=""){
  
  # dfnew <- df[,c("X.N0s_marg","X.N0s_cond")]
  # dfnew <- 100*dfnew
  
  marg <- ggplot(df,aes(100*X.N0s_marg/nrow(df)))+
    geom_histogram(bins=50)+
    scale_x_continuous(limits=c(20,50))+
    ylab("No of populations")+xlab("Percent empty")+
    theme_classic() +
    ggtitle(paste("Marginal -",caption))
  
  cond <- ggplot(df,aes(100*X.N0s_cond/nrow(df)))+
    geom_histogram(bins=50)+
    scale_x_continuous(limits=c(20,50))+
    ylab("No of populations")+xlab("Percent empty")+
    theme_classic() +
    ggtitle(paste("Conditional -",caption))
  
  return(list(marg=marg,cond=cond))
  
}


zerohist_pois_gam00_m0 <- make_zerohist(res_pois_gam00_m0,"gamma=0")
zerohist_pois_gam10_m0 <- make_zerohist(res_pois_gam10_m0,"gamma=-0.1")
zerohist_pois_gam25_m0 <- make_zerohist(res_pois_gam25_m0,"gamma=-0.25")

zerohist_pois_slopes_gam00_m0 <- make_zerohist(res_pois_slopes_gam00_m0,"gamma=0")
zerohist_pois_slopes_gam10_m0 <- make_zerohist(res_pois_slopes_gam10_m0,"gamma=-0.1")
zerohist_pois_slopes_gam25_m0 <- make_zerohist(res_pois_slopes_gam25_m0,"gamma=-0.25")

zerohist_ZIP_gam00_m0 <- make_zerohist(res_ZIP_gam00_m0,"gamma=0")
zerohist_ZIP_gam10_m0 <- make_zerohist(res_ZIP_gam10_m0,"gamma=-0.1")
zerohist_ZIP_gam25_m0 <- make_zerohist(res_ZIP_gam25_m0,"gamma=-0.25")

zerohist_ZIP_slopes_gam00_m0 <- make_zerohist(res_ZIP_slopes_gam00_m0,"gamma=0")
zerohist_ZIP_slopes_gam10_m0 <- make_zerohist(res_ZIP_slopes_gam10_m0,"gamma=-0.1")
zerohist_ZIP_slopes_gam25_m0 <- make_zerohist(res_ZIP_slopes_gam25_m0,"gamma=-0.25")























