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
#' @param complete_effect logical; TRUE performs a (conservative) test whether
#' the whole spline has a significant influence after accounting for all other effects.
#' @details
#' Function provides the test vectors for every location of the given covariate
#'
#'
testvec_for_gamm4 <- function(mod, name, sigma2 = NULL, nrlocs=7,
                              complete_effect = FALSE)
{
  
  if (grepl(",",name)) name <- unlist(strsplit(name, ","))
  modmat <- model.matrix(mod$gam)
  names <- attr(modmat, "dimnames")[[2]]
  xnames <- attr(mod$gam$terms, "term.labels")
  pnames <- c(name, setdiff(xnames, name))
  # ind <- grep(name, names)
  
  if("gamm" %in% class(mod)) {
    sterms <- c(unlist(sapply(mod$gam$smooth,"[[","vn")), names(mod$lme$coefficients$random))
  } else {
    sterms <- c(unlist(sapply(mod$gam$smooth,"[[","vn")), names(mod$mer@cnms))
  }
  SigmaInv <- mod$gam$Vp
  if(is.null(sigma2)) sigma2 <- mod$gam$sig2
  
  if(all(name %in% sterms))
  {
    if(complete_effect){
      # Adapted from iboost package (https://github.com/davidruegamer/iboost/)
      
      # from the selectiveInference package
      svdu_thresh <- function(x) {
        svdx <- svd(x)
        inds <- svdx$d > svdx$d[1] * sqrt(.Machine$double.eps)
        return(svdx$u[, inds, drop = FALSE])
      }
      
      if (length(name) == 1) {
      ind_effect <- which(grepl(paste0("s(",name), colnames(modmat), fixed = TRUE))
      }
      if (length(name) == 2) {
        ind_effect <- which(grepl(paste0("te(",paste(name, collapse = ",")), colnames(modmat), fixed = TRUE))
      }
      
      Xj <- modmat[,ind_effect]
      Xminusj <- svdu_thresh(modmat[,setdiff(1:ncol(modmat),ind_effect)])
      vTs <- list(t(svdu_thresh(Xj - tcrossprod(Xminusj) %*% Xj)))
      attr(vTs[[1]], "var") <- sigma2
      
    }else{ # single point on a non-linear function
      
      if(length(nrlocs)>1) qs <- nrlocs else qs <-
          quantile(mod$gam$model[[name]], seq(0,1,l=nrlocs+2)[-c(1,nrlocs+2)])
      datap <- as.data.frame(cbind(qs, matrix(0, nrow = length(qs),
                                              ncol = length(pnames) - 1)))
      names(datap) <- pnames[apply(sapply(pnames, grepl, x = names), 2, any)]
      lpm <- predict(mod$gam, newdata = datap, type = "lpmatrix")
      # set columns to zero not associated with the variable
      lpm[, which(!grepl(paste0("s(",name), colnames(lpm), fixed = TRUE))] <- 0
      vTs <- lapply(1:nrow(lpm), function(j){
        k <- lpm[j,]%*%SigmaInv%*%t(modmat)/mod$gam$sig2
        attr(k, "var") <- as.numeric(lpm[j,]%*%SigmaInv%*%lpm[j,])*sigma2/mod$gam$sig2
        attr(k, "loc") <- qs[j]
        return(k)
      })
      
    }
    
    #
    # return(lapply(vTs, t))
    
  }else{
    
    # # selind <- which(colnames(modmat)==name)
    # # find name in colnames of model matrix
    # selind <- grepl(pattern = paste0(name, "[:digit:]*"), colnames(modmat)) &
    #   # remove interactions
    #   !grepl(pattern = ":", colnames(modmat), fixed=TRUE)
    # # vTs <- list(((SigmaInv%*%t(modmat))/sigma2)[selind, , drop=FALSE])
    # vTs <- testvec_for_mm(mod = mod$mer, which = which(selind),
    #                       marginal = FALSE, sig2 = sigma2, efficient = TRUE)
    # attr(vTs[[1]], "var") <- (SigmaInv)[selind,selind]*sigma2/mod$gam$sig2
    datap <- as.data.frame(cbind(1, matrix(0, nrow = 1,
                                           ncol = length(pnames) - 1)))
    names(datap) <- pnames[apply(sapply(pnames, grepl, x = names), 2, any)]
    lpm <- predict(mod$gam, newdata = datap, type = "lpmatrix")
    # set columns to zero not associated with the variable
    lpm[, which(!grepl(name, colnames(lpm), fixed = TRUE))] <- 0
    vTs <- lapply(1:nrow(lpm), function(j){
      k <- lpm[j,]%*%SigmaInv%*%t(modmat)/mod$gam$sig2
      attr(k, "var") <- as.numeric(lpm[j,]%*%SigmaInv%*%lpm[j,])*sigma2/mod$gam$sig2
      # attr(k, "loc") <- qs[j]
      return(k)
    })
    
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
#' @param G true random effect covariance
#' @param sig2 true error variance
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
    if(is.null(VCOV)) VCOV <- Matrix::crossprod(C)/sig2 + A
    
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
    vTs <- V0inv%*%Matrix::t(C)/sig2
  }
  
  if(is.null(which)) if(marginal) which <- 1:ncol(X) else which <- 1:ncol(C)
  vTs <- lapply(which, function(j){
    
    vt <- vTs[j,]
    if(!marginal) attr(vt, "var") <- V0inv[j,j]
    return(vt)
    
  })
  return(lapply(vTs,t))
  
}
