#' Function which computes selective p-values and intervals for \code{gamm4} and \code{merMod} objects
#'
#' @param mod an object of class \code{merMod} or result of \code{gamm4} function
#' @param checkFun a function of \code{y}, a vector of the same length as 
#' the original response vector which returns \code{TRUE} or \code{FALSE} depending
#' on whether the selection event for a given \code{y} corresponds to the 
#' original model selection. See the example for more details.
#' @param nrSamples integer; the number of Monte Carlo samples to be used for inference (defaults to 1000)
#' @param bayesian logical; whether or not to use a bayesian type covariance
#' @param sigma2 variance used for inference; per default the estimated variance of \code{mod} is used. 
#' Other options are a conservative estimate based on the variance of the response is used ("varY")
#' or to supply a numeric value to base inference on a customize variance
#' @param VCOV covariance matrix of dimension of the response used for inference; 
#' per default the estimated covariance of \code{mod} is used. 
#' Otherwise a matrix must be supplied on which basis inference is conducted. If the true 
#' covariance is unknown, an conservative alternative to plugging in the estimator is given 
#' by using the covariance of the refitted mixed model, for which all fixed effects but the intercept 
#' are excluded.
#' @param conditional logical; determines whether to use the conditional or marginal approach
#' when \code{mod} is of class \code{merMod}, i.e., inference is sought for a linear mixed model
#' @param name character; for the \code{gamm4}-case: the name of the covariate, for which inference is done
#' @param nrlocs integer; for the \code{gamm4}-case: the number of locations tested for non-linear effects
#' @param which integer; for the \code{merMod}-case: defining the effect for which inference is done
#' @param vT list of vectors (optional); if inference is sought for a customized test vector, this argument
#' can be used
#' @param G true random effect covariance (optional)
#' 
#' @details Note that the additive and conditional mixed model approach currently only works for 
#' a diagonal error covariance.
#' 
#' @examples
#' 
#' library(lme4)
#' library(lmerTest)
#' 
#' ##### BASED ON lmerTest HELP PAGE #########
#' # define function to fit a model based on response
#' modFun <- function(y)
#' {
#'   ham$y <- y
#'   lmer(y ~ Gender + Information * Product + (1 | Consumer) +
#'   (1 | Product), data=ham)
#'   
#'   }
#' 
#' # define a function to select a model (based on a model)
#' selFun <- function(mod) step(mod)
#' 
#' # define a function which extracts the results of the selection procedure
#' extractSelFun <- function(this_mod){
#' 
#' this_mod <- attr(this_mod, "model")
#' if(class(this_mod)=="lm") 
#'   return(attr(this_mod$coefficients, "names")) else
#'     return(c(names(fixef(this_mod)), names(getME(this_mod, "theta"))))
#' 
#' }
#'
#' 
#' ## backward elimination of non-significant effects:
#' (step_result <- selFun(modFun(ham$Informed.liking)))
#' attr(step_result, "model")
#' ## Elimination tables for random- and fixed-effect terms:
#' (sel <- extractSelFun(step_result))
#' 
#' ## Now we can finally define the function checking the congruency
#' ## with the original selection
#' 
#' checkFun <- function(yb){ 
#' 
#'  this_mod <- modFun(yb)
#'  setequal( extractSelFun(selFun(this_mod)), sel )
#'  
#'  }
#' 
#' # Now let's compute valid p-values conditional on the selection
#' res <- mocasin(attr(step_result, "model"), this_y = ham$Informed.liking,
#'           checkFun = checkFun, which = 1:4, nrSamples = 50)
#' 
#' print(res)
#'
mocasin <- function(
  mod, 
  checkFun,
  this_y = NULL,
  nrSamples = 1000,
  bayesian = TRUE,
  VCOV = "est",
  sigma2 = "est",
  conditional  = TRUE,
  name = NULL, 
  nrlocs = 7,
  which = NULL,
  vT = NULL,
  G = NULL,
  efficient = TRUE
)
{
  

  if(is.null(VCOV) | is.null(sigma2)) 
    stop("VCOV and sigma2 must be supplied (see ?mocasin)")
  
  # save variance and covariance if given
  if(!is.character(sigma2)){ 
    this_sigma2 <- sigma2
    sigma2 <- "supplied"
  }
  
  this_VCOV <- NULL
  if(!is.character(VCOV)){ 
    this_VCOV <- VCOV
    VCOV <- "supplied"
  }
  
  # check type
  isMM <- inherits(x = mod, what = "merMod")
  isLM <- class(mod)=="lm"
  
  # type of VCOV
  diagCOV <- (conditional & isMM) | isLM | !isMM
  
  # set bayesian to FALSE for marginal case
  if(!conditional & isMM) bayesian <- FALSE
  
  # get response
  if(is.null(this_y)){
    if(isMM) this_y <- mod@resp else if(isLM)
      this_y <- mod$y else
        this_y <- mod$mer@resp
  }  
  
  # define #obs
  n <- length(this_y)
  
  # define variance
  if(sigma2 == "est")
  {
    
    if(isLM | isMM) sigma2 <- sigma(mod)^2 else
        sigma2 <- mod$gam$sig2
    
  }else if(sigma2 == "varY"){
    
    # plugin by Tibshirani et al. 2018
    sigma2 <- sqrt(var(this_y)*(n-1)/n)    
    
  }else if(sigma2=="supplied"){
    
    sigma2 <- this_sigma2
    
  }else stop("Please supply a numeric value of one of the two options 'varY', 'est' for sigma2.")
  
  # define covariance
  if(VCOV == "est"){
    
    if(diagCOV) # cond. MM, LM or AM
      VCOV <- sigma2 * diag(rep(1,n)) else # marginal MM with plugin
        VCOV <- extract_SigPlZGZT(mod, sig2 = sigma2)
      
  }else if(VCOV == "supplied"){
    
    VCOV <- this_VCOV
    
  }else stop("Please supply a matrix or the option 'est' for VCOV.")
  
  
  if(!is.null(which)) wn <- which else wn <- name
  
  if(isMM){
    
    if(conditional){
      
      vT <- testvec_for_mm(mod, 
                           marginal = !conditional, 
                           G = G, 
                           sig2 = sigma2, 
                           which = wn)

    }else{
      
      # use VCOV only if supplied, else more efficient matrix inversion
      # in the testvec function
      vT <- testvec_for_mm(mod, 
                           marginal = !conditional, 
                           VCOV = this_VCOV,
                           sig2 = sigma2, 
                           which = wn,
                           efficient = efficient)
      
      
      
    }
    
  }else{
    
    
    
    vT <- testvec_for_gamm4(mod, 
                            name = wn, 
                            sigma2 = sigma2, 
                            nrlocs = nrlocs)
    
  }
  
  # compute selective inference
  selinf <- lapply(vT, function(vt)
    pval_vT_cov(vT = vt,
                VCOV = VCOV,
                this_y = this_y,
                nrSamples = nrSamples,
                checkFun = checkFun,
                bayesian = bayesian
    )
  )
  names(selinf) <- wn
  
  retl <- list(vT = vT, selinf = selinf)
  class(retl) <- "selfmade"
  return(retl)
  
}

#' @title Generic methods for selfmade objects
#' 
#' @description Generic methods which can be used for objects fitted with the \code{mocasin} function
#' 
#' @param x selfmade object
#'
#' @method print selfmade
#' @rdname methodsSelfmade
#' 
print.selfmade <- function(x, ...)
{
  
  covariate <- names(x$selinf)
  df <- do.call("rbind", x$selinf)
  df$covariate <- covariate
  print(df)
  invisible(df)
  
}
  