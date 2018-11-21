###########################
##     Report Results    ##
###########################
## 11/21/2018

require(xtable)

## set number of digits for final table
dig <- 2

## load data
## Local:
ests <- data.frame(read.table(file="../Data Analysis/ests.txt",header=T,row.names=1))
SEs <- data.frame(read.table(file="../Data Analysis/SEs.txt",header=T,row.names=1))


## create intervals
CIs_lo <- ests-1.96*SEs
CIs_hi <- ests+1.96*SEs

## exponentiate
ests[1:(nrow(ests)-4),] <- exp(ests[1:(nrow(ests)-4),])
CIs_lo[1:(nrow(CIs_lo)-4),] <- exp(CIs_lo[1:(nrow(CIs_lo)-4),])
CIs_hi[1:(nrow(CIs_hi)-4),] <- exp(CIs_hi[1:(nrow(CIs_hi)-4),])

## make full results table
tab_res <- c()
for(cc in 1:ncol(ests)){ ## loop over analyses
  tab_res <- cbind(tab_res,
                   round(ests[,cc],dig), ## estimates
                   paste0("(",round(CIs_lo[,cc],dig),",",round(CIs_hi[,cc],dig),")") ## CIs
                   )
}

colnames(tab_res) <- paste0(rep(colnames(ests),each=2),
                            rep(c("_Est","_CI"),times=ncol(ests)) )
rownames(tab_res) <- c("e0","e1 DES",
                       "a0","a1 DES","a2 msmk","a3 yob5155","a4 yob5660","a5 yob61plus",
                       "b0","b1 DES","b2 msmk","b3 yob5155","b4 yob5660","b5 yob61plus",
                       "Sigma0","Sigma1",
                       "Gamma0","Gamma1" )


## NAs
tab_res[tab_res=="(NA,NA)"] <- ""
tab_res[is.na(tab_res)] <- ""


# ## all results
# xtable(tab_res)
# ## all marginal results
# xtable(tab_res[,1:(9*2)])
# ## all conditional results
# xtable(tab_res[,(9*2+1):ncol(tab_res)])

## marginal intercepts results
xtable(tab_res[,c(1:10,15:16)])
write.csv(tab_res[,c(1:10,15:16)],"../Data Analysis/marg_int_results.csv")
## marginal slopes results
xtable(tab_res[,c(1:6,11:14,17:18)])
write.csv(tab_res[,c(1:6,11:14,17:18)],"../Data Analysis/marg_slopes_results.csv")

## conditional intercepts results
xtable(tab_res[,18+c(1:2,5:8,13:14)])
write.csv(tab_res[,18+c(1:2,5:8,13:14)],"../Data Analysis/cond_int_results.csv")
## conditional slopes results
xtable(tab_res[,18+c(3:4,9:12,15:16)])
write.csv(tab_res[,18+c(3:4,9:12,15:16)],"../Data Analysis/cond_slopes_results.csv")

