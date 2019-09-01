#' Generate samples based on the decomposition of Y
#'
#' @param orthdir nx1 vector; component orthogonal to the direction of interest
#' @param dir nx1 vector; component in the direction of interest
#' @param this_sd scalar; (estimated / surrogate) value for the root of error variance
#' @param nrSample integer; number of samples to be used
#' @param sampFun function; function to generate samples from
#' @param checkFun function; function to congruency with initial selection
#' @param trace logical; if \code{TRUE}, a progress bar will be printed to the console
#'
gen_samples <- function(
  orthdir,
  dir, 
  this_sd, # DR: why do we need this?
  nrSample = 1000,
  sampFun = function(n) rnorm(n, mean = 0, sd = this_sd),
  checkFun,
  trace = 1)
{
  if(trace)
    pb <- txtProgressBar(min = 0, max = nrSample, style = 3)
  fac <- sampFun(nrSample)
  yb <- lapply(fac, function(tau) as.numeric(orthdir + tau*dir))
  logvals <- c()
  for(i in 1:length(yb)){
    
    logvals[i] <- checkFun(yb[[i]])
    if(trace){
      setTxtProgressBar(pb, i)
    }
    
  }
  if(trace) close(pb)
  
  return(list(logvals = logvals,
              fac = fac))
  
}