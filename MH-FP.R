library(brglm2)

#' Mantel-Haenszel estimate of odds ratio across K 2x2 tables with Firth Penalization
#' @param ctK 2x2xK contingency table
#' @param b0 null value of log OR, default to 0
#' @param alpha significance level, default to 0.05
#' @return a list with: bhat (log OR estimator); se (standard error);  p.val (penalized LRT test and Wald test); ci (confidence interval)
MH.FP <- function(ctK, b0=0, alpha=0.05){
  K = dim(ctK)[3]
  X = rep(1:0, K)
  S = rep(1:K, each=2)
  N1 = as.vector( ctK[,1,] )
  N0 = as.vector( ctK[,2,] )
  N = N0+N1
  ##
  dat = data.frame(N1=N1, N0=N0, X=X, S=factor(S))
  fit1 = glm(cbind(N1,N0)~X+S,family=binomial(logit),method="brglmFit",data=dat)
  X1 = model.matrix(fit1)
  bhat = coef(fit1)
  ## accurate computation of log(1 + e^x)
  L1e = function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  lfn = function(bta) {
    eta = drop(X1 %*% bta)
    p  = plogis(eta)
    ll = sum(N1 * eta - N * L1e(eta))
    jp = 0.5*determinant(crossprod(X1, X1*(N*p*(1-p))), logarithm=TRUE)$modulus
    ans = ll+jp
    return(ans)
  }
  lmax = lfn(bhat)
  ## Inner maximization over nuisance parameters at fixed beta_1
  Pfn = function(b1, start = bhat[-2]){
    negf = function(xpar) { 
      rcf = numeric(ncol(X1)); 
      rcf[2] = b1; 
      rcf[-2] = xpar; 
      -lfn(rcf) 
    }
    ans = optim(start,negf,method="BFGS",control=list(reltol=1e-12,maxit=1000))
    -ans$value
  }
  ### Penalized LR test of H0: beta_1 = b0
  LR = 2*lmax - 2*Pfn(b0)
  Lval = pchisq(LR, df = 1, lower.tail = FALSE)

  ## (1-alpha) profile penalized likelihood CI
  se = sqrt(diag(vcov(fit1)))[2]
  c0 = qchisq(1-alpha, 1)
  fn = function(b1) 2*lmax - 2*Pfn(b1) - c0
  lo = uniroot(fn, c(bhat[2] - 4 * se, bhat[2]), extendInt = "downX")$root
  hi = uniroot(fn, c(bhat[2], bhat[2] + 4 * se), extendInt = "upX")$root
  ci = c(lower = lo, upper = hi)
  pvals = c('PLR'=Lval, 'WT'=2*pnorm(-abs(bhat[2]-b0)/se))
  return(list(bta=bhat[2], se=se, p.val=pvals, ci=ci) )
}
