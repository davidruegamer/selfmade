#' Function which computes selective p-values and intervals for \code{gamm4} and
#' \code{merMod} objects
#'
#' @param mod an object of class \code{merMod} or result of \code{gamm4} function
#' @param checkFun a function of \code{y}, a vector of the same length as
#' the original response vector which returns \code{TRUE} or \code{FALSE} depending
#' on whether the selection event for a given \code{y} corresponds to the
#' original model selection. See the example for more details.
#' @param nrSamples integer; the number of Monte Carlo samples to be used for
#' inference (defaults to 1000)
#' @param bayesian logical; whether or not to use a bayesian type covariance
#' @param conditional logical; determines whether to use the conditional or
#' marginal approach
#' when \code{mod} is of class \code{merMod}, i.e., inference is sought for a
#' linear mixed model
#' @param name character; for the \code{gamm4}-case: the name of the covariate,
#' for which inference is done
#' @param nrlocs integer; for the \code{gamm4}-case: the number of locations
#' tested for non-linear effects
#' @param complete_effect list of logical values for each \code{name}; 
#' TRUE performs a (conservative) test whether
#' the whole spline has a significant influence after accounting for all other effects.
#' @param which integer; for the \code{merMod}-case: defining the effect for
#' which inference is done
#' @param vT list of vectors (optional); if inference is sought for a customized
#' test vector, this argument
#' can be used
#' @param G true random effect covariance (optional)
#' @param trace logical; if TRUE, a progress bar is printed in the console
#' @param this_y original response vector (explicit reference may be necessary
#' for certain model classes)
#' @param varInTestvec for expert use only; variance used in the test vector definition
#' @param varForSampling variance used for inference; per default the estimated variance
#' of \code{mod} is used.
#' Other options are a conservative estimate based on the variance of the
#' response is used ("varY") or to supply a numeric value to base inference
#' on a customize variance
#' @param VCOV_vT for expert use only; VCOV used in the test vector definition
#' @param VCOV_sampling covariance matrix of dimension of the response used for inference;
#' per default the estimated covariance of \code{mod} is used.
#' Otherwise a matrix must be supplied on which basis inference is conducted.
#' If the true
#' covariance is unknown, an conservative alternative to plugging in the
#' estimator is given
#' by using the covariance of the refitted mixed model, for which all fixed
#' effects but the intercept
#' are excluded.
#' @param efficient logical; whether or not to compute the test statistic based on
#' an (efficient) weighted LS estimator instead of a OLS estimator for the marginal model
#'
#'
#' @details Note that the additive and conditional mixed model approach
#' currently only works for a diagonal error covariance.
#'
#' @export
#' @import parallel Matrix lme4
#' @importFrom stats dnorm formula model.matrix predict
#' quantile reformulate rnorm sigma uniroot var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @rawNamespace
#' if(getRversion() >= "3.3.0") {
#' importFrom("stats", sigma)
#' } else {
#' importFrom("lme4", sigma)
#' }
#' @examples
#'
#' library(lme4)
#' if(require(lmerTest)){
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
#'           checkFun = checkFun, which = 1:4, nrSamples = 50, trace = FALSE)
#'
#' # print(res)
#' }
#'
#' # gamm4 example similar to the one from gamm4 help page
#' if(require(gamm4)){
#' set.seed(0)
#' dat <- gamSim(1,n=500,scale=2) ## simulate 4 term additive truth
#'
#' dat$y <- 3 + dat$x0^2 + rnorm(n=500)
#' br <- gamm4(y~ s(x0) + s(x1), data = dat)
#' summary(br$gam) ## summary of gam
#'
#' # do not use any selection
#' checkFun <- function(yb) TRUE
#'
#' res <- mocasin(br, this_y = dat$y,
#'           checkFun = checkFun,
#'           nrlocs = c(0.7,1),
#'           nrSamples = 100)
#' }
mocasin <- function(
  mod,
  checkFun,
  this_y = NULL,
  nrSamples = 1000,
  bayesian = TRUE,
  varInTestvec = c("est", "minMod", "varY", "supplied"),
  varForSampling = c("est", "minMod", "varY", "supplied"),
  VCOV_vT = NULL,
  VCOV_sampling = NULL,
  conditional  = TRUE,
  name = NULL,
  nrlocs = 7,
  complete_effect = NULL,
  which = NULL,
  vT = NULL,
  G = NULL,
  efficient = TRUE,
  trace = TRUE
)
{

  ######### get further setting specitivities #########
  # check type
  isMM <- inherits(x = mod, what = "merMod")
  isLM <- any(class(mod)=="lm")
  isGAMM <- any(class(mod)=="gamm")

  # get response
  if(is.null(this_y)){
    if(isMM) this_y <- mod@resp else if(isLM)
      this_y <- mod$y else
        if(isGAMM) this_y <- mod$lme$data$y else
        this_y <- mod$mer@resp$y
  }

  # match arguments for variance
  varInTestvec = match.arg(varInTestvec)
  varForSampling = match.arg(varForSampling)

  # type of VCOV
  diagCOV <- (conditional & isMM) | isLM | !isMM

  # set bayesian to FALSE for marginal case
  if(!conditional & isMM) bayesian <- FALSE
  #####################################################

  ######### define further models and params ##########
  if(any("minMod" %in% c(varInTestvec, varForSampling))){

    if(!isLM){ # estimate IC mod
      if(isMM) modIC = minMod(mod) else
        modIC = minMod(mod$mer)
    }

    if(isLM){

      sigmaIC2 <- sigma(mod)^2
      tauIC2 <- NULL

    }else{


      sigmaIC2 = sigma(modIC)^2
      tauIC2 = getME(modIC, "theta")

    }

  }


  # define #obs
  n <- length(this_y)

  # estimated variance
  if(isLM | isMM) sigma2 <- sigma(mod)^2 else
    sigma2 <- mod$gam$sig2

  # plugin by Tibshirani et al. 2018
  sigma2_y <- sqrt(var(this_y)*(n-1)/n)

  #####################################################

  ############## define (co-)variances ################

  if(varInTestvec=="supplied" &
     is.null(VCOV_vT))
      stop("Must supply VCOV_vT if varInTestvec=='supplied'")


  if(varForSampling=="supplied" &
    is.null(VCOV_sampling))
      stop("Must supply VCOV_sampling if varForSampling=='supplied'")

  # if no supplied variance for testvector, define

  sigma2_sampling <- NULL

  if(is.null(VCOV_sampling)){

    if(conditional){

      VCOV_sampling <- switch (varForSampling,
                               est = sigma2,
                               varY = sigma2_y,
                               minMod = sigmaIC2
      )

      VCOV_sampling = VCOV_sampling * diag(rep(1,n))

    }else{

      if(varForSampling == "varY")
        stop("Option 'varY' not possible for marginal mixed model inference.")

      VCOV_sampling <- switch(varForSampling,
                              est = mod@sigma^2 * diag(rep(1,n)),
                              minMod = vcov_RI(modIC)
      )

      sigma2_sampling <- switch(varForSampling,
                                varY = sigma2_y,
                                est = mod@sigma^2,
                                minMod = modIC@sigma^2)

    }
  }

  if(is.null(sigma2_sampling)) sigma2_sampling = VCOV_sampling[1,1]

  # if no supplied variance for testvector, define
  if(is.null(VCOV_vT)){

    if(conditional){

      VCOV_vT <- switch (varInTestvec,
                         est = NULL,
                         varY = sigma2_y,
                         minMod = sigmaIC2
      )

    }else{ # marginal



    }

  }

  #####################################################

  #### define testvector corresponding to settings ####

  if(!is.null(which)) wn <- which else wn <- name

  if(isMM){

    if(conditional){

      vT <- testvec_for_mm(mod,
                           marginal = !conditional,
                           G = G,
                           sig2 = VCOV_vT,
                           which = wn)

    }else{ # marginal


      # use VCOV only if supplied, else more efficient matrix inversion
      # in the testvec function
      vT <- testvec_for_mm(mod,
                           marginal = !conditional,
                           VCOV = VCOV_vT,
                           sig2 = sigma2_sampling,
                           which = wn,
                           efficient = efficient)



    }

  }else{ # additive case


    if(is.null(wn)) wn <- attr(mod$gam$terms, "term.labels")
    if(is.null(complete_effect)){ 
      complete_effect <- rep(FALSE, length(wn))
      names(complete_effect) <- wn
    }
    
    vT <- sapply(wn, function(name)
      res <- testvec_for_gamm4(mod,
                               name = name,
                               sigma2 = VCOV_vT,
                               nrlocs = nrlocs,
                               complete_effect[[name]])
    )
    
    if(is.list(vT[[1]])) vT <- unlist(vT, recursive = FALSE)
    
  }
  #####################################################

  ############## compute selective inference ##########
  pbi <- 0
  selinf <- lapply(1:length(vT), function(j){

    vt = vT[[j]]
    ce = complete_effect[[wn[j]]]
    if(trace) cat("Computing inference for variable (location) ", pbi+1, "\n\n")
    res <- pval_vT_cov(vT = vt,
                       VCOV = VCOV_sampling,
                       this_y = this_y,
                       nrSamples = nrSamples,
                       checkFun = checkFun,
                       bayesian = bayesian,
                       trace = trace,
                       complete_effect = ce
    )
    if(trace) cat("\n\n")
    pbi <<- pbi + 1
    return(res)
  }
  )
  # if(!isMM) names(selinf) <- wn
  
  ######################################################

  retl <- list(vT = vT, selinf = selinf)
  class(retl) <- "selfmade"
  return(retl)

}

#' @title Generic methods for selfmade objects
#'
#' @description Generic methods which can be used for objects
#' fitted with the \code{mocasin} function
#'
#' @param x selfmade object
#' @param ... further arguments, currently unused.
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
