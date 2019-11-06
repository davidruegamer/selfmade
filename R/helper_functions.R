
minMod <- function(mod)
{
  
  lhs <- formula(mod)[[2]]
  reFormula <- cnms2formula(mod@cnms)
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

# copied from cAIC4 package
cnms2formula <-
  function(cnms) {
    # A function that builds a random effects formula from the ?component names?,
    # a list that can be extracted from an lmerMod object by .@cnms or
    # getME(., "cnms").
    #
    # Args:
    #   cnms     = List from an lmerMod object by .@cnms or getME(., "cnms").
    #
    # Returns:
    #   reFormula = random effects part of a lmerMod formula
    #
    len      <- unlist(lapply(cnms, length))
    cnms     <- cnms[which(len != 0)]
    charForm <- character(length(cnms))
    
    for(i in 1:length(cnms)) {
      if (cnms[[i]][1] == "(Intercept)") {
        cnms[[i]][1] <- "1"
      } else {
        tpv <- cnms[[i]]
        cnms[[i]] <- append("",tpv)
        cnms[[i]][1] <- "-1"
      }
      charForm[i] <- paste("(", paste(cnms[[i]], collapse = " + "),
                           " | ",names(cnms)[i], ")", sep = "")
    }
    
    reFormula <- paste(charForm, collapse = " + ")
    
    return(reFormula)
  }


hatmatfun_gamm <- function(obj, 
                           nr_smooths=length(obj$gam$smooth), 
                           what=c("Vmat", "hatmat", "coefmat")
){
  
  # special hat matrix for gamms
  # see https://researchportal.bath.ac.uk/files/9228764/tgamm4.pdf
  # page 15
  
  what <- match.arg(what)
  
  lmeobj <- obj$lme
  gamobj <- obj$gam
  y <- gamobj$y
  X <- predict(gamobj, type = "lpmatrix")
  sigma2 <- gamobj$sig2
  
  Vlme <- extract.lme.cov2(lmeobj,
                           gamobj$model,
                           nr_smooths+1)
  Vlme_ind <- Vlme$ind
  Vlme <- Vlme$V
  if(is.list(Vlme)) Vlme <- Matrix::bdiag(Vlme)
  
  nr_fixef <- length(attr(gamobj$pterms,"term.labels")) + 
    attr(mod$gam$terms,"intercept")
  P <- Matrix::bdiag(append(rep(0,nr_fixef),
                            lapply(1:length(gamobj$smooth),
                                   function(i) gamobj$sp[i] * 
                                     Reduce(sum, gamobj$smooth[[i]]$S))))
  bayV <- gamobj$Vp
  
  # this is not very stable. Should switch to the approach by mgcv
  # in the near future.
  vlmeinv = solve(Vlme)
  Vmat <- solve(t(X[Vlme_ind,])%*%vlmeinv%*%X[Vlme_ind,] + P/sigma2)
  if(what=="Vmat") return(Vmat)
  coefmat <- Vmat%*%t(X)%*%vlmeinv
  if(what=="coefmat") return(coefmat)
  if(what=="hatmat") return(X%*%coefmat) else return(coefmat%*%y)
  
}