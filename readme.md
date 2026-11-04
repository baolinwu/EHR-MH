# MH-FP
 - Estimating Risk Factor Odds Ratios for Rare Events from Stratified EHR Summary Data
 - Reference
   - Wu et. al (2026) Assessing risk factors for rare clinical events from summary-level electronic health records: a confounder-adjusted approach. tech rep.
 - R scripts for conducting analyses and simulation studies
```R
library(brglm2)
library(logistf)
source("MH-FP.R")    # provides MH.FP()

## Simulated proctectomy summary data: one row per stratum defined by
## surgeon type (GS = 1 for general surgeon) and the covariates;
## Y1 = number of deaths, Y0 = number of survivors
dat1 = read.table("MortalityCount.txt", header = TRUE)
dat1$XMale = 1 - dat1$XFemale
fml = ~ GS + Age1 + Age2 + RaceBlack + RaceWhite + XMale

## LR: maximum likelihood, likelihood ratio inference
lr = glm(update(fml, cbind(Y1, Y0) ~ .), family = binomial, data = dat1)
drop1(lr, test = "LRT")        ## likelihood ratio test for each term
confint(lr, "GS")              ## profile likelihood CI

## LR-FP: Firth's penalized likelihood
flr = glm(update(fml, cbind(Y1, Y0) ~ .), family = binomial, method = "brglmFit", data = dat1)
summary(flr)                   ## estimates (Wald inference)
## profile penalized likelihood inference, using frequency weights
dat2 = rbind(transform(dat1, y = 1, w = Y1), transform(dat1, y = 0, w = Y0))
dat2 = dat2[dat2$w > 0, ]
flr.pl = logistf(update(fml, y ~ .), data = dat2, weights = w)
summary(flr.pl)

## 2 x 2 x 30 table: rows = (GS, CRS), columns = (death, survival),
## layers = covariate strata
strat = interaction(dat1$Age1, dat1$RaceBlack, dat1$RaceWhite, dat1$XMale, drop = TRUE)
tbi = array(0, dim = c(2, 2, nlevels(strat)))
for(k in seq_len(nlevels(strat))) {
  s = dat1[strat == levels(strat)[k], ]
  tbi[1, , k] = c(sum(s$Y1[s$GS == 1]), sum(s$Y0[s$GS == 1]))
  tbi[2, , k]= c(sum(s$Y1[s$GS == 0]), sum(s$Y0[s$GS == 0]))
}

## MH (CMH test, default with continuity correction) and MH-exact
mantelhaen.test(tbi)
mantelhaen.test(tbi, exact = TRUE)

## MH-FP, with all strata and excluding strata without deaths
MH.FP(tbi)
MH.FP(tbi, b0=0.5) ## test H_0: \beta=0.5
MH.FP(tbi[, , apply(tbi[, 1, ], 2, sum) > 0])
 

## load the simulation scenario 4 data: 
  ## X: main exposure of interest; covariates: Xage, Xage2, RaceBlack, RaceWhite, Xsex
  ## S: covariate strata id (1-30)
  ## NN: strata specific total counts (sample size)
  ## prob0, prob1, prob2: strata specific event probabilities (corresponding to log OR: 0, 0.4, 0.8)
load('data_sc4.rda')
## load(url('https://github.com/baolinwu/EHR-MH/blob/main/data_sc4.rda?raw=true'))
## load other scenarios: load('data_sc1.rda'); load('data_sc2.rda'); load('data_sc3.rda')
 K = dim(data_sc4)[1]
 ## simulate events under null of OR=1
 Y = rbinom(K, data_sc4$NN,data_sc4$prb0)
 ## simulate events under log OR = 0.4
 Y = rbinom(K, data_sc4$NN,data_sc4$prb1)
 ## simulate events under log OR = 0.8
 Y = rbinom(K, data_sc4$NN,data_sc4$prb2)
 dati = data.frame(data_sc4,Y1=Y,Y0=data_sc4$NN-Y)
 ##
 flr = glm(cbind(Y1,Y0) ~ X+Xage+Xage2+RaceBlack+RaceWhite+Xsex, family=binomial(logit), method="brglmFit", data=dati)
 summary(flr)
 ##
 dati2 = rbind(transform(dati, y = 1, w = Y1), transform(dati, y = 0, w = Y0))
 dati2 = dati2[dati2$w > 0, ]
 flr.pl = logistf(y ~ X+Xage+Xage2+RaceBlack+RaceWhite+Xsex, data = dati2, weights = w)
 summary(flr.pl)
## 2x2x30 contingency table
 tbi = array(0, dim=c(2,2,30))
 for(k in 1:30){
    tbi[2:1,1,k] = dati$Y1[2*(k-1)+1:2]
    tbi[2:1,2,k] = dati$Y0[2*(k-1)+1:2]
 }    
 mantelhaen.test(tbi)
 mantelhaen.test(tbi, exact=TRUE)
 MH.FP(tbi)

## check OR
 data_sc4$prb2[2*(1:30)]/(1-data_sc4$prb2[2*(1:30)])/(data_sc4$prb2[2*(1:30)-1]/(1-data_sc4$prb2[2*(1:30)-1])) ## OR
 data_sc4$prb2[2*(1:30)]/data_sc4$prb2[2*(1:30)-1] ## RR ~ OR
```
