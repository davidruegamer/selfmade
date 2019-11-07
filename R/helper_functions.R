
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
                           what=c("Vmat", "hatmat", "coefmat","coef")
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
  
  # Switch to the more efficient version?
  Vlme <- extract.lme.cov(lmeobj,
                          gamobj$model,
                          nr_smooths+1)
  
  nr_fixef <- length(attr(gamobj$pterms,"term.labels")) + 
    attr(mod$gam$terms,"intercept")
  # P <- Matrix::bdiag(append(rep(0,nr_fixef),
  #                           lapply(1:length(gamobj$smooth),
  #                                  function(i) gamobj$sp[i] * 
  #                                    Reduce(sum, gamobj$smooth[[i]]$S))))
  
  first.para <- sapply(1:length(gamobj$smooth), function(i) gamobj$smooth[[i]]$first.para)
  last.para <- sapply(1:length(gamobj$smooth), function(i) gamobj$smooth[[i]]$last.para)
  nxf <- ncol(X)
  k <- 1
  
  S <- matrix(0, nxf, nxf)
  
  for (i in 1:length(gamobj$smooth)) {
    if (is.null(gamobj$smooth[[i]]$fac)) {
      ind <- first.para[i]:last.para[i]
      ns <- length(gamobj$smooth[[i]]$S)
      if (ns) 
        for (l in 1:ns) {
          S[ind, ind] <- S[ind, ind] + gamobj$smooth[[i]]$S[[l]] * 
            gamobj$sp[k]
          k <- k + 1
        }
    }
    else {
      flev <- levels(gamobj$smooth[[i]]$fac)
      ind <- first.para[i]:(first.para[i] + gamobj$smooth[[i]]$n.para - 
                              1)
      ns <- length(gamobj$smooth[[i]]$S)
      for (j in 1:length(flev)) {
        if (ns) 
          for (l in 1:ns) {
            S[ind, ind] <- S[ind, ind] + gamobj$smooth[[i]]$S[[l]] * 
              gamobj$sp[k]
            k <- k + 1
          }
        k <- k - ns
        ind <- ind + gamobj$smooth[[i]]$n.para
      }
      k <- k + ns
    }
  }
  
  
  bayV <- gamobj$Vp
  
  # this is not very stable. Should switch to the approach by mgcv
  # in the near future.
  vlmeinv = solve(Vlme)
  Vmat <- solve(t(X)%*%vlmeinv%*%X + S/sigma2)
  if(what=="Vmat") return(Vmat)
  coefmat <- Vmat%*%t(X)%*%vlmeinv
  if(what=="coefmat") return(coefmat)
  if(what=="hatmat") return(X%*%coefmat) else return(coefmat%*%y)
  
}
