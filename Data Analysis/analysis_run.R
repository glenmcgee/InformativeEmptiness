#######################
##    NHS Analysis   ##
#######################
## 11/21/2018
## Fit models to NHS data
## prepare dataframes for JMMICS pack
## then fit:
 ## Conditional:
  ## joint conditional (no 0s)
  ## joint conditional
  ## outcome-only glmm
 ## Marginal:
  ## GEE-Exch
  ## IEE
  ## WGEE
  ## JMM
  ## marginal size model (for comparing joint models)

## load packages
source("./Functions/JMMICS.R")
## load and clean data
source("./Data Analysis/analysis_clean_data.R")
## prelim tables
source("./Data Analysis/analysis_EDA.R")
## other
require(tidyr)
require(dplyr)
require(glmmML)
require(lme4)
require(geepack)
require(ggplot2)
require(xtable)
library(JMMICSpack)



################################
### prepare data for JMMICSpack
YY <- datG2$adhd
XX <- cbind(1,datG2$desqx1,datG2$msmk2,datG2$yob89_5155,datG2$yob89_5660,datG2$yob89_61plus)#,datG2$raceWhite,datG2$momed2,datG2$momed3,datG2$momed4)
NN <- datG1$totalkids
ZZ <- cbind(1,datG1$desqx1,datG1$msmk2,datG1$yob89_5155,datG1$yob89_5660,datG1$yob89_61plus)
ZZ0 <- ZZ[,1:2]
IDD <- datG2$id2
### data for JMMICSpack (no 0s)
YY_no0 <- datG2_no0$adhd
XX_no0 <- cbind(1,datG2_no0$desqx1,datG2_no0$msmk2,datG2_no0$yob89_5155,datG2_no0$yob89_5660,datG2_no0$yob89_61plus)
NN_no0 <- datG1_no0$totalkids
ZZ_no0 <- cbind(1,datG1_no0$desqx1,datG1_no0$msmk2,datG1_no0$yob89_5155,datG1_no0$yob89_5660,datG1_no0$yob89_61plus)
IDD_no0 <- datG2_no0$id2




########################################
##            Model Fits              ##
########################################

########################################
## joint (marginal) models
##
gee <- geeglm(adhd~desqx1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2,family=binomial,corstr="exchangeable")
gee_est <- gee$coef
gee_SE <- summary(gee)$coef[,2]
print("gee")
##
iee <- geeglm(adhd~desqx1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2,family=binomial,corstr="independence")
iee_est <- iee$coef
iee_SE <- summary(iee)$coef[,2]
print("iee")
##
wgee <- geeglm(adhd~desqx1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2,family=binomial,corstr="independence",weights=(1/totalkids))
wgee_est <- wgee$coef
wgee_SE <- summary(wgee)$coef[,2]
print("wgee")
##
jmm_no0s <- JMMICS_fit(Nk=NN_no0,Zk=ZZ_no0,Yki=YY_no0,Xki=XX_no0,IDk=IDD_no0,Z0k=ZZ_no0,weights=NA,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad=50)
jmm_no0s_est <- jmm_no0s[[1]]
jmm_no0s_SE <- jmm_no0s[[2]]
print("jmm_no0s")
##
jmm <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad=50)
jmm_est <- jmm[[1]]
jmm_SE <- jmm[[2]]
print("jmm")
##
jmm_slopes_no0s <- JMMICS_fit(Nk=NN_no0,Zk=ZZ_no0,Yki=YY_no0,Xki=XX_no0,IDk=IDD_no0,Z0k=ZZ_no0,weights=NA,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad=50)
jmm_slopes_no0s_est <- jmm_slopes_no0s[[1]]
jmm_slopes_no0s_SE <- jmm_slopes_no0s[[2]]
print("jmm_slopes_no0s")
##
jmm_slopes <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad=50)
jmm_slopes_est <- jmm_slopes[[1]]
jmm_slopes_SE <- jmm_slopes[[2]]
print("jmm_slopes")
##
jmm_ZIP <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=TRUE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad=50)
jmm_ZIP_est <- jmm_ZIP[[1]]
jmm_ZIP_SE <- jmm_ZIP[[2]]
print("jmm_ZIP")
##
jmm_ZIP_slopes <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=TRUE,slopes=TRUE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,nquad=50)
jmm_ZIP_slopes_est <- jmm_ZIP_slopes[[1]]
jmm_ZIP_slopes_SE <- jmm_ZIP_slopes[[2]]
print("jmm_ZIP_slopes")

########################################
## joint (conditional) models
##
naive <- glmmML(adhd~desqx1+msmk2+yob89_5155+yob89_5660+yob89_61plus,cluster=id2,family=binomial,data=datG2,method="ghq",n.points=30)
naive_est <- c(naive$coefficients,naive$sigma)
naive_SE <- c(naive$coef.sd,naive$sigma.sd)
print("naive")
##
naive_slopes <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=FALSE,nquad=50)
naive_slopes_est <- naive_slopes[[1]]
naive_slopes_SE <- naive_slopes[[2]]
print("naive_slopes")
##
joint_no0s <- JMMICS_fit(Nk=NN_no0,Zk=ZZ_no0,Yki=YY_no0,Xki=XX_no0,IDk=IDD_no0,Z0k=ZZ_no0,weights=NA,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad=50)
joint_no0s_est <- joint_no0s[[1]]
joint_no0s_SE <- joint_no0s[[2]]
print("joint_no0s")
##
joint <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad=50)
joint_est <- joint[[1]]
joint_SE <- joint[[2]]
print("joint")
##
joint_slopes_no0s <- JMMICS_fit(Nk=NN_no0,Zk=ZZ_no0,Yki=YY_no0,Xki=XX_no0,IDk=IDD_no0,Z0k=ZZ_no0,weights=NA,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad=50)
joint_slopes_no0s_est <- joint_slopes_no0s[[1]]
joint_slopes_no0s_SE <- joint_slopes_no0s[[2]]
print("joint_slopes_no0s")
##
joint_slopes <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad=50)
joint_slopes_est <- joint_slopes[[1]]
joint_slopes_SE <- joint_slopes[[2]]
print("joint_slopes")
##
joint_ZIP <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=TRUE,slopes=FALSE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad=50)
joint_ZIP_est <- joint_ZIP[[1]]
joint_ZIP_SE <- joint_ZIP[[2]]
print("joint_ZIP")
##
joint_ZIP_slopes <- JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ0,weights=NA,minNk=0,NegBin=FALSE,ZIP=TRUE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=TRUE,joint=TRUE,nquad=50)
joint_ZIP_slopes_est <- joint_ZIP_slopes[[1]]
joint_ZIP_slopes_SE <- joint_ZIP_slopes[[2]]
print("joint_ZIP_slopes")


########################################
##          Collect Results           ##
########################################
n_a <- ncol(ZZ)        ## no. of alphas
n_b <- ncol(XX)        ## no. of betas
n_e <- ncol(ZZ0)       ## no. of epsilons

## order of parameters:  epsilon,       alphas,        betas,    sigma0, sigma1,   gamma0, gamma1 
                    ## c(rep(NA,n_e), rep(NA,n_a),  rep(NA,n_b),   rep(NA,2),        rep(NA,2))
  ## where sigma0=sigma for random intercepts

## order that parameters are output from JMMICSpack: alpha gamma0 gamma1 beta sigma0 sigma1 epsilon

#c(rep(NA,n_e),rep(NA,n_a),rep(NA,n_b),rep(NA,2),rep(NA,2)),

ests <- cbind(c(rep(NA,n_e),rep(NA,n_a),gee_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_e),rep(NA,n_a),iee_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_e),rep(NA,n_a),wgee_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_e),jmm_no0s_est[1:n_a],jmm_no0s_est[(n_a+1)+(1:n_b)],jmm_no0s_est[n_a+1+n_b+1],NA,jmm_no0s_est[(n_a+1)],NA),
              c(rep(NA,n_e),jmm_est[1:n_a],jmm_est[(n_a+1)+(1:n_b)],jmm_est[(n_a+1+n_b)+1],NA,jmm_est[n_a+1],NA),
              c(rep(NA,n_e),jmm_slopes_no0s_est[1:n_a],jmm_slopes_no0s_est[(n_a+2)+(1:n_b)],jmm_slopes_no0s_est[(n_a+2+n_b)+(1:2)],jmm_slopes_no0s_est[(n_a)+(1:2)]),
              c(rep(NA,n_e),jmm_slopes_est[1:n_a],jmm_slopes_est[(n_a+2)+(1:n_b)],jmm_slopes_est[(n_a+2+n_b)+(1:2)],jmm_slopes_est[(n_a)+(1:2)]),
              c(jmm_ZIP_est[(n_a+1+n_b+1)+1:n_e],jmm_ZIP_est[1:n_a],jmm_ZIP_est[(n_a+1)+(1:n_b)],jmm_ZIP_est[(n_a+1+n_b)+1],NA,jmm_ZIP_est[n_a+1],NA),
              c(jmm_ZIP_slopes_est[(n_a+2+n_b+2)+1:n_e],jmm_ZIP_slopes_est[1:n_a],jmm_ZIP_slopes_est[(n_a+2)+(1:n_b)],jmm_ZIP_slopes_est[(n_a+2+n_b)+(1:2)],jmm_ZIP_slopes_est[(n_a)+(1:2)]),
              c(rep(NA,n_e),rep(NA,n_a),naive_est[1:n_b],naive_est[(n_b)+1],NA,rep(NA,2)),
              c(rep(NA,n_e),rep(NA,n_a),naive_slopes_est[1:n_b],naive_slopes_est[(n_b)+(1:2)],rep(NA,2)),
              c(rep(NA,n_e),joint_no0s_est[1:n_a],joint_no0s_est[(n_a+1)+(1:n_b)],joint_no0s_est[n_a+1+n_b+1],NA,joint_no0s_est[(n_a+1)],NA),
              c(rep(NA,n_e),joint_est[1:n_a],joint_est[(n_a+1)+(1:n_b)],joint_est[(n_a+1+n_b)+1],NA,joint_est[n_a+1],NA),
              c(rep(NA,n_e),joint_slopes_no0s_est[1:n_a],joint_slopes_no0s_est[(n_a+2)+(1:n_b)],joint_slopes_no0s_est[(n_a+2+n_b)+(1:2)],joint_slopes_no0s_est[(n_a)+(1:2)]),
              c(rep(NA,n_e),joint_slopes_est[1:n_a],joint_slopes_est[(n_a+2)+(1:n_b)],joint_slopes_est[(n_a+2+n_b)+(1:2)],joint_slopes_est[(n_a)+(1:2)]),
              c(joint_ZIP_est[(n_a+1+n_b+1)+1:n_e],joint_ZIP_est[1:n_a],joint_ZIP_est[(n_a+1)+(1:n_b)],joint_ZIP_est[(n_a+1+n_b)+1],NA,joint_ZIP_est[n_a+1],NA),
              c(joint_ZIP_slopes_est[(n_a+2+n_b+2)+1:n_e],joint_ZIP_slopes_est[1:n_a],joint_ZIP_slopes_est[(n_a+2)+(1:n_b)],joint_ZIP_slopes_est[(n_a+2+n_b)+(1:2)],joint_ZIP_slopes_est[(n_a)+(1:2)])
              )
              
SEs <- cbind(c(rep(NA,n_e),rep(NA,n_a),gee_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_e),rep(NA,n_a),iee_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_e),rep(NA,n_a),wgee_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_e),jmm_no0s_SE[1:n_a],jmm_no0s_SE[(n_a+1)+(1:n_b)],jmm_no0s_SE[n_a+1+n_b+1],NA,jmm_no0s_SE[(n_a+1)],NA),
              c(rep(NA,n_e),jmm_SE[1:n_a],jmm_SE[(n_a+1)+(1:n_b)],jmm_SE[(n_a+1+n_b)+1],NA,jmm_SE[n_a+1],NA),
              c(rep(NA,n_e),jmm_slopes_no0s_SE[1:n_a],jmm_slopes_no0s_SE[(n_a+2)+(1:n_b)],jmm_slopes_no0s_SE[(n_a+2+n_b)+(1:2)],jmm_slopes_no0s_SE[(n_a)+(1:2)]),
              c(rep(NA,n_e),jmm_slopes_SE[1:n_a],jmm_slopes_SE[(n_a+2)+(1:n_b)],jmm_slopes_SE[(n_a+2+n_b)+(1:2)],jmm_slopes_SE[(n_a)+(1:2)]),
              c(jmm_ZIP_SE[(n_a+1+n_b+1)+1:n_e],jmm_ZIP_SE[1:n_a],jmm_ZIP_SE[(n_a+1)+(1:n_b)],jmm_ZIP_SE[(n_a+1+n_b)+1],NA,jmm_ZIP_SE[n_a+1],NA),
              c(jmm_ZIP_slopes_SE[(n_a+2+n_b+2)+1:n_e],jmm_ZIP_slopes_SE[1:n_a],jmm_ZIP_slopes_SE[(n_a+2)+(1:n_b)],jmm_ZIP_slopes_SE[(n_a+2+n_b)+(1:2)],jmm_ZIP_slopes_SE[(n_a)+(1:2)]),
              c(rep(NA,n_e),rep(NA,n_a),naive_SE[1:n_b],naive_SE[(n_b)+1],NA,rep(NA,2)),
              c(rep(NA,n_e),rep(NA,n_a),naive_slopes_SE[1:n_b],naive_slopes_SE[(n_b)+(1:2)],rep(NA,2)),
              c(rep(NA,n_e),joint_no0s_SE[1:n_a],joint_no0s_SE[(n_a+1)+(1:n_b)],joint_no0s_SE[n_a+1+n_b+1],NA,joint_no0s_SE[(n_a+1)],NA),
              c(rep(NA,n_e),joint_SE[1:n_a],joint_SE[(n_a+1)+(1:n_b)],joint_SE[(n_a+1+n_b)+1],NA,joint_SE[n_a+1],NA),
              c(rep(NA,n_e),joint_slopes_no0s_SE[1:n_a],joint_slopes_no0s_SE[(n_a+2)+(1:n_b)],joint_slopes_no0s_SE[(n_a+2+n_b)+(1:2)],joint_slopes_no0s_SE[(n_a)+(1:2)]),
              c(rep(NA,n_e),joint_slopes_SE[1:n_a],joint_slopes_SE[(n_a+2)+(1:n_b)],joint_slopes_SE[(n_a+2+n_b)+(1:2)],joint_slopes_SE[(n_a)+(1:2)]),
              c(joint_ZIP_SE[(n_a+1+n_b+1)+1:n_e],joint_ZIP_SE[1:n_a],joint_ZIP_SE[(n_a+1)+(1:n_b)],joint_ZIP_SE[(n_a+1+n_b)+1],NA,joint_ZIP_SE[n_a+1],NA),
              c(joint_ZIP_slopes_SE[(n_a+2+n_b+2)+1:n_e],joint_ZIP_slopes_SE[1:n_a],joint_ZIP_slopes_SE[(n_a+2)+(1:n_b)],joint_ZIP_slopes_SE[(n_a+2+n_b)+(1:2)],joint_ZIP_slopes_SE[(n_a)+(1:2)])
             )


rownames(ests) <- rownames(SEs) <- c("e0","e1 DES",
                                     "a0","a1 DES","a2 msmk","a3 yob5155","a4 yob5660","a5 yob61plus",
                                     "b0","b1 DES","b2 msmk","b3 yob5155","b4 yob5660","b5 yob61plus",
                                     "Sigma0","Sigma1",
                                     "Gamma0","Gamma1"
                                     )
colnames(ests) <- colnames(SEs) <- c("GEE","IEE","WEE",
                                     "JMM No0s","JMM",
                                     "JMMSlopes No0s","JMMSlopes",
                                     "JMM ZIP","JMMSlopes ZIP",
                                     "Out-Only","Out-Only Slopes",
                                     "Joint No0s","Joint",
                                     "JointSlopes No0s","JointSlopes",
                                     "Joint ZIP","JointSlopes ZIP" 
                                     )


write.table(ests,file="~/Data Analysis/ests.txt")
write.table(SEs,file="~/Data Analysis/SEs.txt")






xtable(ests)
xtable(SEs)



