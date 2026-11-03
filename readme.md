# MH-FP
 - Estimation of odds ratio across K 2x2 tables with sparse counts
 - Sample R codes
```R
 library(brglm2)
 ## analyze the "sample proctectomy" data
 dat1 = read.table('MortalityCount.txt', header=TRUE)
 dat1$XMale = 1- dat1$XFemale
 lr = glm(cbind(Y1,Y0) ~ GS+Age1+Age2+RaceBlack+RaceWhite+XMale, family=binomial(logit), data=dat1)
 summary(lr)
 flr = glm(cbind(Y1,Y0) ~ GS+Age1+Age2+RaceBlack+RaceWhite+XMale, family=binomial(logit), method="brglmFit", data=dat1)
 summary(flr)
 ## 2x2x30 contingency table
 tbi = array(0, dim=c(2,2,30))
 for(k in 1:30){
    tbi[2:1,1,k] = dat1$Y1[2*(k-1)+1:2]
    tbi[2:1,2,k] = dat1$Y0[2*(k-1)+1:2]
 }    
 mantelhaen.test(tbi)
 mantelhaen.test(tbi, exact=TRUE)
 mantelhaen.test(tbi, correct=FALSE)
 MH.FP(tbi)
 ##
 j1 = which(tbi[1,1,] + tbi[2,1,]>0);  length(j1)
 MH.FP(tbi[,,j1])
 
 ## load the simulation scenario 4 data: 
  ## X is the main exposure of interest; covariates: Xage, Xage2, RaceBlack, RaceWhite, Xsex
 load('data_sc4.rda')
 K = dim(data_sc4)[1]
 ## simulate events under null of OR=1
 Y = rbinom(K, data_sc4$NN,data_sc4$prb0)
 ## simulate events under log OR = 0.4
 Y = rbinom(K, data_sc4$NN,data_sc4$prb1)
 ## simulate events under log OR = 0.8
 Y = rbinom(K, data_sc4$NN,data_sc4$prb2)
 dati = data.frame(data_sc4,Y1=Y,Y0=data_sc4$NN-Y)
 flr = glm(cbind(Y1,Y0) ~ X+Xage+Xage2+RaceBlack+RaceWhite+Xsex, family=binomial(logit), method="brglmFit", data=dati)
 summary(flr)
 ## 2x2x30 contingency table
 tbi = array(0, dim=c(2,2,30))
 for(k in 1:30){
    tbi[2:1,1,k] = dati$Y1[2*(k-1)+1:2]
    tbi[2:1,2,k] = dati$Y0[2*(k-1)+1:2]
 }    
 mantelhaen.test(tbi)
 mantelhaen.test(tbi, exact=TRUE)
 mantelhaen.test(tbi, correct=FALSE)
 MH.FP(tbi)
```
