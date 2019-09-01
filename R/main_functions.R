#' Calculate inference based on survived samples and weights
#'
#' @param survr vector of survived samples
#' @param tstat the original test statistic
#' @param w weights for the survived samples
#' @param var_est (estimated) variance of the test statistic
#' @param alpha see \code{\link{mocasin}}
#'
selinf <- function(
  survr,
  tstat,
  w,
  var_est,
  alpha
)
{

  # check if any sample has survived
  nrsurv <- length(survr)
  if(nrsurv==0) return(data.frame(nrsurv = nrsurv, pval = 0, cil = -Inf, ciu = Inf))

  # which sample are more extreme than the observed value
  l2 <- survr > tstat
  survr_gr <- survr[l2]
  survr_le <- survr[!l2]

  # calculate the one- and two-sided p-value
  # if(all(w==0)) pv <- 0 else{
    pe <- sum(w[l2]) / sum(w)
    pv <- 2*min(pe, 1-pe)
  # }

  # calculate the CIs
  ftlo <- function(t) sum(exp(survr_gr * t / (var_est[1])) * w[l2]) /
    sum(exp(survr * t / (var_est[2])) * w)
  ftup <- function(t) sum(exp(survr_le * t / (var_est[1])) * w[!l2]) /
    sum(exp(survr * t / (var_est[2])) * w)

  testvals <- seq(min(survr) - 8*sqrt(var_est[1]),
                  max(survr) + 8*sqrt(var_est[1]),
                  l = 1000)
  flovals <- sapply(testvals, ftlo)
  fupvals <- sapply(testvals, ftup)
  ll <- min(which(!is.na(flovals) &
                    !is.nan(flovals)))
  lu <- max(which(!is.na(flovals) &
                    !is.nan(flovals)))

  low <- try(uniroot(f = function(x) ftlo(x) - alpha/2,
                     interval = testvals[c(ll,lu)],
                     extendInt = "upX")$root)

  up <- try(uniroot(f = function(x) ftup(x) - alpha/2,
                    interval = testvals[c(ll,lu)],
                    extendInt = "downX")$root)

  if(class(low)=="try-error") low <- -Inf
  if(class(up)=="try-error") up <- Inf

  ci <- c(low, up)

  return(data.frame(nrsurv = nrsurv, tstat = tstat, pval = pv, cil = low, ciu = up))

}


#' Calculate selective p-value for given covariance
#'
#' @param vT test vector of function
#' @param VCOV covariance used for distribution of test statistic
#' @param this_y original response vector
#' @param nrSamples number of samples to be used
#' @param checkFun function; function to congruency with initial selection
#' @param twosided logical; compute two-sided p-values?
#' @param bayesian see \code{\link{mocasin}}
#' @param alpha see \code{\link{mocasin}}
#' @param maxiter maximum number of iteratoins to perform the linesearch used
#' in the sampling procedure
#' @param trace see \code{\link{mocasin}}
#' @param complete_effect logical; TRUE performs a (conservative) test whether
#' the whole spline has a significant influence after accounting for all other effects.
#' @param ... further arguments passed to vT if specified as function
#'
pval_vT_cov <- function(
  vT,
  VCOV,
  this_y,
  nrSamples,
  checkFun,
  twosided = TRUE,
  bayesian = FALSE,
  alpha = 0.05,
  maxiter = 10,
  trace = TRUE,
  complete_effect = FALSE,
  ...
)
{
  
  if(complete_effect){
    # group test
    # Adapted from iboost package (https://github.com/davidruegamer/iboost/)
    
    m <- nrow(vT)
    
    if(m==1) vT <- vT / as.numeric(sqrt(tcrossprod(vT)))
    R <- vT %*% this_y
    Z <- sqrt(sum(R^2))
    u <- t(vT) %*% R / Z
    yperp <- this_y - u * Z
    var <- attr(vT, "var")
    
    sigma1 <- sqrt(m - 2 * (gamma((m+1)/2) / gamma(m/2))^2 )
    
    # Do Monte Carlo (Importance Sampling)
    ss <- gen_samples(checkFun = checkFun, 
                      this_sd = NULL,
                      sampFun = function(B){ 
                        
                        rBs <- Z + rnorm(B) * sigma1 * sqrt(var)
                        rBs[rBs>0]
                        
                      }, 
                      nrSample = nrSamples, 
                      orthdir = yperp, 
                      dir = u)
    
    rBs <- ss$fac
    survr <- rBs[ss$logvals]
    
    Z <- Z/sqrt(var)
    s <- survr/sqrt(var)
    w <- (s^(m-1)) * (exp(-s^2/2)) / dnorm(s, mean = Z, sd = sigma1)
    var_est <- var
    tstat <- Z
    
    
    
  }else{ # standard "univariate" test
    
    n <- length(this_y)
    vvT <- tcrossprod(vT)
    tstat <- as.numeric(vT%*%this_y)
    
    if(!bayesian) var_est <- as.numeric(vT%*%VCOV%*%t(vT)) else
      var_est <- attr(vT, "var")
    dirV <- (VCOV%*%t(vT)/var_est)
    orthdir <- (diag(n) - dirV%*%vT)%*%this_y
    
    samples <- gen_samples(
      orthdir = orthdir,
      dir = dirV,
      this_sd = sqrt(var_est),
      sampFun = function(n) rnorm(n, mean = tstat, sd = sqrt(var_est)),
      nrSample = nrSamples,
      checkFun = checkFun,
      trace = trace)
    
    # extract survived samples and weights
    survr <- samples$fac[samples$logvals]
    nom <- dnorm(survr, mean = 0, sd = sqrt(var_est))
    denom <- dnorm(survr, mean = tstat, sd = sqrt(var_est))
    
    var_est <- rep(var_est, 2)
    
    while(sum(nom)==0 & all(denom!=0) & maxiter-1 > 0){
      
      var_est[2] <- var_est[2] * abs(tstat)/sqrt(var_est[2])
      
      samples <- gen_samples(
        orthdir = orthdir,
        dir = dirV,
        this_sd = sqrt(var_est[1]),
        sampFun = function(n) rnorm(n, mean = tstat, sd = var_est[2]),
        nrSample = nrSamples,
        checkFun = checkFun)
      
      survr <- samples$fac[samples$logvals]
      nom <- dnorm(survr, mean = 0, sd = sqrt(var_est[1]))
      denom <- dnorm(survr, mean = tstat, sd = var_est[2])
      
      maxiter <- maxiter - 1
      
    }
    
    w <- nom / denom

  }

  # compute p-value and CI
  return(
    selinf(
      survr = survr,
      tstat = tstat,
      w = w,
      var_est = var_est,
      alpha = alpha
    )
  )


}


