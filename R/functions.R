#' Generate samples based on the decomposition of Y
#' 
#' @param orthdir nx1 vector; component orthogonal to the direction of interest
#' @param dir nx1 vector; component in the direction of interest
#' @param this_sd scalar; (estimated / surrogate) value for the root of error variance
#' @param nrSample integer; number of samples to be used
#' @param sampFun function; function to generate samples from
#' @param checkFun function; function to congruency with initial selection
#' 
gen_samples <- function(
  orthdir, 
  dir,
  this_sd,
  nrSample = 1000,
  sampFun = function(n) rnorm(n, mean = 0, sd = this_sd),
  checkFun)
{
  
  fac <- sampFun(nrSample)
  yb <- sapply(fac, function(tau) as.numeric(orthdir + tau*dir))
  logvals <- apply(yb, 2, checkFun)
  
  return(list(logvals = logvals,
              fac = fac))
  
}


#' Calculate inference based on survived samples and weights
#' 
#' @param survr vector of survived samples
#' @param tstat the original test statistic
#' @param w weights for the survived samples
#' @param var_est (estimated) variance of the test statistic
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
  pe <- sum(w[l2]) / sum(w)
  pv <- 2*min(pe, 1-pe)
  
  # calculate the CIs
  ftlo <- function(t) sum(exp(survr_gr * t / (var_est)) * w[l2]) / 
    sum(exp(survr * t / (var_est)) * w)
  ftup <- function(t) sum(exp(survr_le * t / (var_est)) * w[!l2]) / 
    sum(exp(survr * t / (var_est)) * w)
  
  testvals <- seq(min(survr) - 8*sqrt(var_est), 
                  max(survr) + 8*sqrt(var_est),
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
  
  return(data.frame(nrsurv = nrsurv, pval = pv, cil = low, ciu = up))
  
}


#' Calculate selective p-value for given covariance
#' 
#' @param vT test vector of function
#' @param VCOV covariance used for distribution of test statistic
#' @param this_y original response vector
#' @param coefNr which coefficient p-values should be calculated for
#' @param nrSamples number of samples to be used
#' @param checkFun function; function to congruency with initial selection
#' @param twosided logical; compute two-sided p-values?
#' @param bayesian_sd if not NULL, this variance is used as variance of the test statistic.
#' The idea of this argument is to allow for usage of a bayesian covariance, but 
#' this argument can also be used otherwise.
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
  ...
)
{
  
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
    checkFun = checkFun)
  
  # extract survived samples and weights
  survr <- samples$fac[samples$logvals]
  w <- dnorm(survr, mean = 0, sd = sqrt(var_est)) / 
    dnorm(survr, mean = tstat, sd = sqrt(var_est))
  
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




#' Function to calculate the test vector for an object fitted with \code{gamm4}
#'
#' @param mod an object fitted with \code{gamm4}
#' @param name character; name of the covariate for which inference should be
#' calculated
#' @param sigma2 variance to be used in the covariance definition. If \code{NULL},
#' the estimate \code{mod$gam$sig2} is used.
#' @param nrlocs number of locations at which p-values and intervals are to be computed 
#' for non-linear terms. This directly corresponds to a sequence of \code{nrlocs} quantiles 
#' of the given covariate values.
#' @details 
#' Function provides the test vectors for every location of the given covariate
#' 
#'
testvec_for_gamm4 <- function(mod, name, sigma2 = NULL, nrlocs=7) #, plot=FALSE, truefun=NULL, this_y=NULL)
{
  
  modmat <- model.matrix(mod$gam)
  names <- attr(modmat, "dimnames")[[2]]
  xnames <- attr(mod$gam$terms, "term.labels")
  pnames <- c(name, setdiff(xnames, name))
  ind <- grep(name, names)
  sterms <- sapply(mod$gam$smooth,"[[","vn")
  SigmaInv <- mod$gam$Vp
  if(is.null(sigma2)) sigma2 <- mod$gam$sig2
  
  if(name %in% sterms)
  {
    if(length(nrlocs)>1) qs <- nrlocs else qs <- 
        quantile(mod$gam$model[[name]], seq(0,1,l=nrlocs)[-c(1,nrlocs)])
    datap <- as.data.frame(cbind(qs, matrix(0, nrow = length(qs), 
                                            ncol = length(pnames) - 1)))
    names(datap) <- pnames[apply(sapply(pnames, grepl, x = names), 2, any)]
    lpm <- predict(mod$gam, newdata = datap, type = "lpmatrix")
    vTs <- lapply(1:nrow(lpm), function(j){ 
      k <- lpm[j,]%*%SigmaInv%*%t(modmat)/sigma2
      attr(k, "var") <- as.numeric(lpm[j,]%*%SigmaInv%*%lpm[j,])
      attr(k, "loc") <- qs[j]
      return(k)
    })
    # 
    # return(lapply(vTs, t))

  }else{
    
    selind <- which(colnames(modmat)==name)
    vTs <- list((SigmaInv%*%t(modmat)/sigma2)[selind, , drop=FALSE])
    attr(vTs[[1]], "var") <- (SigmaInv)[selind,selind]
    
  }
  
  return(vTs)
  
  # if(plot)
  # {
  #   
  #   plot(mod$gam, select = which(name%in%attr(mod$gam$terms,"term.labels")), ylim = c(-7,9))
  #   intcpt <- coef(mod$gam)[1]
  #   points(qs, sapply(vTs, function(vT) vT%*%this_y) - intcpt)
  #   points(qs, truefun(qs) - intcpt,  col="red")
  #   
  # }

}

#' Function to calculate the test vector for an \code{merMod} object
#'
#' @param mod \code{merMod} object 
#' @param which integer; if NULL, test vector is created for all 
#' coefficients, otherwise for the which'th coefficient
#' @param marginal logical; whether to construct a test vector
#' from the marginal or condition perspective
#' @param VCOV covariance used for the construction of test vector 
#' @param efficient logical; whether to compute the 
#' test vector corresponding to the efficient beta estimator
#' in the marginal case
#'
#' @details 
#' Function provides 2 different marginal test vectors
#' (efficient = TRUE / FALSE) and one conditional test vector
#' 
#' Note that the covariance of residuals R is assumed to be diagonal.
#' 
#' 
#'
#'
testvec_for_mm <- function(
  mod, 
  which = NULL, 
  marginal = TRUE, 
  VCOV = NULL,
  G = NULL,
  sig2 = NULL,
  efficient = TRUE
  )
{
  
  if(class(mod)=="lm"){
    
    X <- model.matrix(mod)
    vTs <- solve(crossprod(X))%*%t(X)
    if(is.null(which)) which <- 1:ncol(X)
    return(lapply(which, function(j) t(vTs[j,])))
    
  } 
  if(is.null(sig2)) sig2 <- sigma(mod)^2
  
  Z       <- getME(mod, "Z")
  X       <- getME(mod, "X")
  n       <- nrow(X)
  Lambda <- getME(mod, "Lambda")
  Lambdat <- getME(mod, "Lambdat")
  
  
  if(!marginal){ 
    
    if(length(mod@Gp)>3) stop("Please implement for more than 2 REs.")
    C <- cbind(X,Z)
    if(is.null(G) || dim(G)[2] != dim(getME(mod,"ST")[[1]])[2]) 
      A <- bdiag(list(matrix(0, ncol=ncol(X), nrow=ncol(X)), Lambda%*%Lambdat)) else
        A <- bdiag(list(matrix(0, ncol=ncol(X), nrow=ncol(X)), 
                        bdiag(list(solve(G))[rep(1,ncol(Z)/ncol(G))])))
    if(class(A)=="try-error") # G does not fit to the selected model
      A <- bdiag(list(matrix(0, ncol=ncol(X), nrow=ncol(X)), Lambda%*%Lambdat))
    if(is.null(VCOV)) VCOV <- crossprod(C)/sig2 + A
    
  }
               
               
  if(marginal){
    
    L       <- getME(mod, "L")
    if(is.null(VCOV)){
      
      V0inv   <- diag(rep(1, n)) - crossprod(solve(L, system = "L") %*% 
                                               solve(L, Lambdat, system = "P") %*% t(Z))
      # = solve(diag(n)+Z%*%Lambda%*%Lambdat%*%t(Z))
      V0inv <- V0inv / sig2

    }else{
      
      V0inv <- solve(VCOV)
      
    }
      
  }  

  if(marginal){
    if(efficient)
      vTs <- solve(crossprod(X, V0inv%*%X))%*%t(X)%*%V0inv else
        vTs <- solve(crossprod(X))%*%t(X)
  }else{ # conditional
    V0inv <- solve(VCOV)
      vTs <- V0inv%*%t(C)/sig2 
  }
  
  if(is.null(which)) if(marginal) which <- 1:ncol(X) else which <- 1:ncol(C)
  vTs <- lapply(which, function(j){ 
    
    vt <- vTs[j,]
    if(!marginal) attr(vt, "var") <- V0inv[j,j]
    return(vt)
  
  })
  return(lapply(vTs,t))
  
}

minMod <- function(mod)
{
  
  lhs <- formula(mod)[[2]]
  reFormula <- cAIC4:::cnms2formula(mod@cnms)
  #modIC <- 
  update(mod, formula. = reformulate(c("1", reFormula), lhs))
  # return(list(sigma = getME(modIC, "sigma"),
  #             tau = getME(modIC, "theta"))
  # )
  
}

vcov_RI <- function(mod, sigma = getME(mod, "sigma"), tau = getME(mod, "theta"),
                    sigmaT = NULL, tauT = NULL)
{
 
  if(!is.null(sigmaT)) sigma <- sigmaT
  if(!is.null(tauT)) tau <- tauT
  n <- NROW(mod@resp$y)
  nrsubj <- nlevels(mod@flist[[1]])
  bdiag(list(diag(n/nrsubj)*sigma^2 + tau^2)[rep(1,nrsubj)]) 
  
}

extract_SigPlZGZT <- function(mod, sig2 = NULL) 
{
 
  if(is.null(sig2)) sig2 <- sigma(mod)^2
  n <- length(mod@response)
  sig2 * (diag(n) + crossprod(getME(mod, "Lambdat") %*% getME(mod, "Zt")))
  
}