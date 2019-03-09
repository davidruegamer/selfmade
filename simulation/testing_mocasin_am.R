# Load libraries and source functions
# TODO: Clean up and use library for selfmade
library(lme4)
library(cAIC4)
library(gamm4)
library(Matrix)
library(MASS)
library(ggplot2)
library(parallel)
library(reshape2)
library(mvtnorm)
source("functions.R")
source("selfmade.R")

set.seed(0) 

n = 500
nrSamples = 500
nrSimIter = 500
saveResults = TRUE

### create data
centeredRN <- function() as.numeric(scale(rnorm(n), scale=F))

x1 <- centeredRN()
x2 <- centeredRN()
x3 <- centeredRN()
x4 <- centeredRN()
x5 <- centeredRN()
# x6 <- rnorm(n)
dat <- data.frame(x1 = x1, x2 = x2, x3 = x3,
                  x4 = x4, x5 = x5) # , x6 = x6)
### define dgp
dgp <- function(n, sd) 3 + dat$x1^2 + sin(dat$x2) + rnorm(n=n, sd)

### define selection functions
critFun = AIC
compareFun = which.min
selfun = function(y) 
{
  
  dat$y <- y
  br <- gam(y ~ x1 + x2 + s(x3) + x4, data = dat)
  br2 <- gam(y ~ s(x1) + x2 + x3 + x4, data = dat)
  br3 <- gam(y ~ s(x1) + s(x2) + x3 + x4, data = dat)
  br4 <- gam(y ~ s(x1) + x2 + s(x3) + x4, data = dat)
  br5 <- gam(y ~ x1 + x2 + x3 + x4, data = dat)
  ret = compareFun(sapply(list(br, br2, br3, br4, br5), critFun))
  attr(ret, "form") = list(br, br2, br3, br4, br5)[[ret]]$formula
  return(ret)
}

###############################################################
################ saving options and filenames #################
###############################################################


saveName = "test_mocasin_gamAIC2_"
regexpr = ".*AIC2\\_(.*)\\.RDS"
plotName = "ggres_test_mocasin2"
plot_final_result = TRUE
nrCores = 50

###############################################################
################# settings for simulation #####################
###############################################################


settings = expand.grid(
  modelType = c("mixed model", "additive model"),
  bayesian = c(TRUE, FALSE),
  conditional = c(TRUE, FALSE),
  efficient = c(TRUE, FALSE),
  varInTestvec = c("est", "minMod", "varY", "supplied"),
  varForSampling = c("est", "minMod", "varY", "supplied"),
  sd = c(1,sqrt(10))
)


### redefine settings (do not use this for overall simulation)
settings$marginal = settings$modelType == "mixed model" & 
  !settings$conditional

drop = 
  (settings$modelType == "additive model" & !settings$conditional) | 
  (!settings$efficient & !settings$marginal) |
  (settings$modelType == "additive model" & (settings$varInTestvec == "minMod" | 
                                               settings$varForSampling == "minMod")) #| 

settings <- settings[!drop,]
# set the model type manually -- CHANGE THIS
settings <- settings[settings$modelType=="additive model",]

# set factors to character
settings[,sapply(settings, class)=="factor"] <- 
  lapply(settings[,sapply(settings, class)=="factor"],
         as.character)
VCOV_vT = VCOV_sampling = NULL


####################################################################
############### Start iterating over simulations ###################
####################################################################


for(i in 1:nrow(settings)){
  
  if(settings[i,"varInTestvec"]=="supplied"){
    # supply true value
    VCOV_vT = settings[i,"sd"]^2
    
  }
  
  if(settings[i,"varForSampling"]=="supplied"){
    # supply true value
    VCOV_sampling = settings[i,"sd"]^2 * diag(rep(1,n))
    
  }
  
  res <- mclapply(1:nrSimIter, function(j){
    
    st <- Sys.time()
    set.seed(j)
    dat$y <- dgp(n=n, sd=settings[i,"sd"])
    # r <- NULL
    # if(checkFun(dat$y))
    wm = selfun(dat$y)
    checkFun <- function(yb) selfun(yb)==wm
    r <- mocasin(mod = gamm4(formula = attr(wm,"form"), data = dat), 
                 this_y = dat$y,
                 checkFun = checkFun, 
                 nrlocs = c(0.7,1), 
                 nrSamples = nrSamples,
                 bayesian = settings[i,"bayesian"],
                 conditional = settings[i,"conditional"],
                 efficient = settings[i,"efficient"],
                 varInTestvec = settings[i,"varInTestvec"],
                 varForSampling = settings[i,"varForSampling"],
                 VCOV_sampling = VCOV_sampling,
                 VCOV_vT = VCOV_vT)   
    print(Sys.time()-st)
    r$sel = wm
    return(r)
    
  }, mc.cores = nrCores)
  
  if(saveResults) saveRDS(res, file=paste0(saveName,i,".RDS"))
  
  res <- res[!sapply(res, is.null)]
  res <- do.call("rbind", lapply(res, function(r){
    
    sel <- r$sel
    r <- r$selinf
    nr <- names(r)
    r <- do.call("rbind", lapply(r, function(x) x[,c("pval"),drop=F]))
    r$var = nr
    r$sel = sel
    rmelt = melt(r, id.vars = c("sel", "var"))[,c(1,2,4)]
    return(rmelt)
    
  }))
  
}

### plot result in the end

if(plot_final_result){
  
  lf <- list.files(pattern = saveName, full.names = T)
  lres <- lapply(lf, function(l){
    
    d <- readRDS(l)
    d <- d[!sapply(d,is.null)]
    sels <- unlist(sapply(d,"[[","sel"))
    d <- do.call("rbind", lapply(1:length(sels), function(i){
      var <- names(d[[i]]$selinf)
      data <- do.call("rbind", d[[i]]$selinf)
      data$sel <- sels[i]
      data$var <- var
      return(data)
    }))
    d$setting_nr = gsub(regexpr,"\\1",l)
    return(d)
    
  })
  
  res <- do.call("rbind", lres)
  
  settings$setting_nr = 1:nrow(settings)
  res <- merge(res, settings, by="setting_nr")
  
  res$sel <- as.factor(res$sel)
  res$varForSampling <- as.factor(res$varForSampling)
  res$varInTestvec <- as.factor(res$varInTestvec)
  
  print(
    gg <- ggplot(data = res,
                 aes(sample = pval,
                     colour = interaction(varInTestvec, varForSampling),
                     linetype = bayesian)) +
      geom_abline(slope = 1, intercept = 0, linetype=2) +
      geom_qq(#,
        distribution = stats::qunif, size = 0.3, geom="line") +
      theme_bw() +# coord_flip() +
      xlab("Expected Quantiles") + ylab("Observed p-values") +
      facet_grid(sel ~ var)
  )
  
  
  saveRDS(gg, file=paste0(plotName,".RDS"))
  
}