##################
##    NHS EDA   ##
##################
## 11/21/2018
## Generate tables and plots for EDA of NHS data
## 1 - generate table of ADHD and DES by cluster size
## 2 - generate table of emptiness rates in each covariate category
## 3 - some prelim tables for ADHD/DES rates by covariates
## 4 - run  preliminary regressions on size, emptiness, and outcome to determine presence of ICS, etc.

## load packages
require(tidyr)
require(dplyr)
require(xtable)
## load clean data
source("../Data Analysis/analysis_clean_data.R")
 

######################
## some summary stats
## no. nurses
length(unique(dat$id))
## rate of adhd diagnosis
mean(datG2$adhd)
## rate of "having at least one child with adhd given that you had at least one child"
mean((datG1_no0$adhd>0))
## rate of multiple diagnoses, given at least one diagnosis
mean((datG1$adhd*datG1$totalkids>1)[datG1_no0$adhd>0])


###################
## 1 Nk table
#### collapse large cluster sizes
dat$Nk_collapse5plus <- dat$totalkids; dat$Nk_collapse5plus[dat$Nk_collapse5plus>=5] <- "5+"
datG1$Nk_collapse5plus <- datG1$totalkids; datG1$Nk_collapse5plus[datG1$Nk_collapse5plus>=5] <- "5+"
tabNk_1 <- dat %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                     Out_pct=round(100*mean(adhd),2))

tabNk_2 <- datG1 %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                     Exp_pct=round(100*mean(desqx1),2))
tabNk <- cbind(tabNk_1,tabNk_2[,3])
tabNk[is.na(tabNk)] <- "--"
colnames(tabNk) <- c("N_k","No. G1s","% ADHD","% DES")
xtable(tabNk,include.rownames = FALSE)
write.table(tabNk,"../Data Analysis/tab1.csv")

################################################
## 2 Table of overall, mothers, and non-mothers
## setting rate=TRUE gives the rates of emptiness in each covariate category ## This one probably more useful
## setting rate=FALSE just gives overall cell percentages
nm_table <- function(EDAvar="N0",cats,names=NA,df=datG1,G1s=TRUE,digits=0,rate=FALSE){
  if(length(names)<length(cats)){
    rnames <- c("Overall",cats)
  }else{
    rnames <- c("Overall",names)
  }
  
  ## column of variable of interest
  EDAvec <- df[,EDAvar]
  
  ## first get overall number of G1s
  total <- c(nrow(df))
  EDA_1 <- c(sum(EDAvec!=1)) ## N>=1
  EDA_0 <- c(sum(EDAvec==1)) ##N=0

  
  ## get category specific numbers
  for(var in cats){
    varvec <- df[,var]
    total <- c(total,sum(varvec==1))
    EDA_1 <- c(EDA_1,sum((EDAvec!=1) & (varvec==1)))
    EDA_0 <- c(EDA_0,sum((EDAvec==1) & (varvec==1)))
  }
  
  
  if(rate==FALSE){
    # percentage
    total_pct <- round(100*total/nrow(df),digits)
    EDA_1_pct <- round(100*EDA_1/nrow(df),digits)
    EDA_0_pct <- round(100*EDA_0/nrow(df),digits)
    ## make table
    tab <- cbind(total,total_pct,EDA_1,EDA_1_pct,EDA_0,EDA_0_pct)
    colnames(tab) <- c("Overall No.","Overall %","N>=1 No.","N>=1 %","N=0 No.","N=0 %")
    
  }else {
    # percentage
    EDA_1_pct <- round(100*EDA_1/sum(EDAvec!=1),digits)
    EDA_0_pct <- round(100*EDA_0/sum(EDAvec==1),digits)  
    ## make table
    tab <- cbind(total,EDA_1,EDA_1_pct,EDA_0,EDA_0_pct)
    colnames(tab) <- c("Overall No.","N>=1 No.","N>=1 Rate %","N=0 No.","N=0 Rate %")
  }
  
  
  #
  rownames(tab) <- rnames
  return(tab)
  
}

nmcats <- c("desqx1","desqx0",
           "msmk2","msmk0",
           "yob89_4650","yob89_5155","yob89_5660","yob89_61plus",
           "raceWhite","raceAA","raceAsian","raceOther",
           "hisp89","nonhisp",
           "momed2","momed3","momed4","momed5","momed6",
           "daded2","daded3","daded4","daded5","daded6",
           "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
           "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
           "phome05","nophome")

nmnames <- c("DES Yes","DES No",
              "Smoke Yes","Smoke No",
              "YOB 45-60","51-55","56-60","61+",
              "White","African American","Asian","Other",
              "Hispanic","Nonhisp",
              "MomEd 2","MomEd 3","MomEd 4 ","MomEd 5 ","MomEd 6 ",
              "DadEd 2","DadEd 3","DadEd 4 ","DadEd 5 ","DadEd 6 ",
              "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
              "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
              "PHOME","PHOME No")

## table of overall/mothers/nonmothers (cell percentages)
nonmothers_tab <- nm_table("N0",cats=nmcats,names=nmnames,df=datG1,digits=1)
xtable(nonmothers_tab,digits=c(0,0,1,0,1,0,1))
## table of overall/mothers/nonmothers (rate of emptiness) ## this one is probably more meaninful if emptiness is the focus
nonmothers_rate_tab <- nm_table("N0",cats=nmcats,names=nmnames,df=datG1,digits=1,rate=TRUE)
xtable(nonmothers_rate_tab,digits=c(0,0,0,1,0,1))



###########################
#### 3 EDA Table  
## (less complete but still interesting plots of DES and ADHD rates by covariate categories)
EDA_table <- function(EDAvar,cats,names=NA,df=datG1,G1s=TRUE,digits=0){
  if(length(names)<length(cats)){
    rnames <- c("Overall",cats)
  }else{
    rnames <- c("Overall",names)
  }
  
  ## column of variable of interest
  EDAvec <- df[,EDAvar]
  
  ## first get overall number of G1s
  EDA_K <- c(sum(EDAvec))
  total_K <- c(nrow(df))
  
  ## get category specific numbers
  for(var in cats){
    total_K <- c(total_K,sum(df[,var]==1))
    EDA_K <- c(EDA_K,sum(EDAvec[df[,var]==1]))
  }
  
  # percentage
  EDA_pct <- round(100*EDA_K/total_K,digits)
  
  ## make table
  tab <- cbind(total_K,EDA_K,EDA_pct)
  rownames(tab) <- rnames
  if(G1s==TRUE){
    colnames(tab) <- c("No. of G1s","No. Exposed","%")
  }else{
    colnames(tab) <- c("No. of G2s","No. Exposed","%")
  }
  
  
  return(tab)
  
}

DEScats<-c("N0","N1","N2","N3","N4plus",
           "yob89_4650","yob89_5155","yob89_5660","yob89_61plus",
           "scand89","ocauc89","afric89","hisp89","asian89","oanc89",
           "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
           "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
           "momed2","momed3","momed4","momed5","momed6",
           "daded2","daded3","daded4","daded5","daded6",
           "msmk2","msmk3","phome05")
DESnames <- c("Nk=0","1","2","3","4+",
              "YOB 45-60","51-55","56-60","61+",
              "Scand","Ocauc","Afric","Hisp","Asian","Oanc",
              "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
              "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
              "MomEd 2","MomEd 3","MomEd 4 ","MomEd 5 ","MomEd 6 ",
              "DadEd 2","DadEd 3","DadEd 4 ","DadEd 5 ","DadEd 6 ",
              "MSMK 2","MSMK 3","PHOME")

## DES table (G1 level)
DES_tab <- EDA_table("desqx1",cats=DEScats,names=DESnames,df=datG1,G1s=TRUE,digits=2)

## ADHD table (G2 level)
ADHDcats <- c("desqx1",DEScats)
ADHDnames <- c("DES",DESnames)
ADHD_tab <- EDA_table("adhd",cats=ADHDcats,names=ADHDnames,df=dat[dat$N0==0,],G1s=FALSE,digits=2)



####################################
## 4 Prelim Regressions

###################################
## function to clean up reg results
clean_reg <- function(g,glmm=FALSE,digits=2){
  if(glmm==FALSE){
    coef <- exp(g$coef)
    # conf <- exp(confint(g))
    SE <- sqrt(diag(summary(g)$cov.scaled))
    conf <- exp(cbind(g$coef-1.96*SE,g$coef+1.96*SE))
  }else{
    coef <- exp(summary(g)$coefficients[,1])
    SE <- summary(g)$coefficients[,2]
    conf <- exp(cbind(summary(g)$coefficients[,1]-1.96*SE,summary(g)$coefficients[,1]+1.96*SE))
  }
  coef <- round(coef,digits)
  conf <- round(conf,digits)
  conf <- paste("(",conf[,1],", ",conf[,2],")",sep="")
  return(cbind(coef,conf))
}


#############################################
## preliminary size regressions
g_pois <- glm(totalkids~desqx1+msmk2+adhd+yob89_5155+yob89_5660+yob89_61plus,data=datG1,family=poisson)
g_nonzero <- glm((1-N0)~desqx1+msmk2+adhd+yob89_5155+yob89_5660+yob89_61plus,data=datG1,family=poisson)

results_Nkreg <- cbind(clean_reg(g_pois),clean_reg(g_nonzero))
colnames(results_Nkreg) <- c("RR","CI","OR","CI")
xtable(results_Nkreg)


