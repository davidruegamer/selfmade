# Load libraries and source functions
# TODO: Clean up and use library for selfmade
library(lme4)
library(cAIC4)
library(lmerTest)
library(gamm4)
library(Matrix)
library(MASS)
library(ggplot2)
library(parallel)
library(reshape2)
library(mvtnorm)
library(selfmade)

set.seed(0) 

n = 500
nrSamples = 500
nrSimIter = 500
saveResults = TRUE


n <- 100
p <- 6
nrsubj <- 5

X <- cbind(rep(1,n), matrix(rnorm(n*(p)), ncol = p))

beta <- c(2, 4, -2, 1, rep(0, p-3))
subjind <- gl(n = nrsubj, k = n / nrsubj)

etaFix <- X%*%beta
signal <- sd(etaFix)

colnames(X) <- c("Intercept", paste("V", 1:6, sep=""))

dat <<- data.frame(X, subjind)

###############################################################
################ saving options and filenames #################
###############################################################


saveName = "test_mocasin_mmcStep_"
regexpr = ".*mmcStep\\_(.*)\\.RDS"
plotName = "ggres_test_mocasin_mmc"
plot_final_result = FALSE
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
  SNR = c(2,4),
  tau = c(4)
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
settings <- settings[settings$modelType!="additive model",]
settings$varInTestvec <- NULL
settings <- settings[!(!settings$conditional & settings$bayesian),]
settings <- settings[!(settings$varForSampling=="varY" & 
                         !settings$conditional),]
settings <- unique(settings)

# set factors to character
settings[,sapply(settings, class)=="factor"] <- 
  lapply(settings[,sapply(settings, class)=="factor"],
         as.character)
VCOV_vT = VCOV_sampling = NULL
settings$sd = signal / settings$SNR

####################################################################
############### Start iterating over simulations ###################
####################################################################


for(i in 1:nrow(settings)){
  
  if(settings[i,"varForSampling"]=="supplied"){
    # supply true value
    VCOV_sampling = settings[i,"sd"]^2 * diag(rep(1,n))
    
  }
  
  res <- mclapply(1:nrSimIter, function(j){
    
    st <- Sys.time()
    sigma <- settings$sd[i]
    tau <- settings$tau[i]
    covTaus <- 0.5*(0.5*tau) # = correlation * sqrt(tau * tauSlope)
    
    vcovTRUE <- matrix(c(tau, covTaus, covTaus, 0.5 * tau), ncol = 2)
    cholVCT <- chol(vcovTRUE)/sigma
    theta <- c(cholVCT)[-2]
    names(theta)<-c("subjind.(Intercept)",
                    "subjind.V2.(Intercept)",
                    "subjind.V2")
    
    form <- as.formula(
      c("~ V1 + V2 + V3 + V4 + V5 + V6 + (1 + V2|subjind)")
    )
    set.seed(j)
    this_y <- simulate(form, newdata = dat, family = gaussian,
                  newparams = list(theta = theta, beta = beta, 
                                   sigma = settings[i,"sd"]))[[1]]
    
    
    ### define selection functions
    modFun <- function(yy) 
    {
      dat$y <- as.numeric(yy)
      lmer(y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + 
             (1 + V2|subjind), REML = FALSE, data = dat)
    }
    
    selFun <- function(mod)
    {
      
      attr(step(mod, reduce.random = FALSE), "model")
      
    }
    
    extractSelFun <- function(this_mod){
      
      if(class(this_mod)=="lm") 
        return(attr(this_mod$coefficients, "names")) else
          return(c(names(fixef(this_mod)), names(getME(this_mod, "theta"))))
      
    }
    
    
    selected_model = selFun(modFun(this_y))
    selection = extractSelFun(selected_model)
    # define function which checks congruency
    checkFun <- function(yb){ 
      
      setequal( extractSelFun(selFun(modFun(yy = yb))), 
                selection )
      
    }
    nameFun = function(this_mod){ 
      
      sel_vars = setdiff(extractSelFun(this_mod), 
                         c("subjind.(Intercept)",
                           "subjind.V2.(Intercept)", 
                           "subjind.V2"))
      if("(Intercept)" %in% sel_vars & 
         "V1" %in% sel_vars & 
         "V2" %in% sel_vars & 
         "V3" %in% sel_vars){
        
        # get inference for V1 if true model selected
        # and else (if upper model), get 
        if(length(sel_vars)==4) return(2) else
          return(5)
        
      }else{
        
        # not correct model selected
        if(length(sel_vars)>1) return(2) else return(1)
        
      }
      
    }
    which = nameFun(selected_model)
    if(is.null(VCOV_sampling)){
      
      
      
    }
    r <- mocasin(mod = selected_model, 
                 this_y = this_y,
                 checkFun = checkFun, 
                 nrSamples = nrSamples,
                 which = which,
                 bayesian = settings[i,"bayesian"],
                 conditional = settings[i,"conditional"],
                 efficient = settings[i,"efficient"],
                 varForSampling = settings[i,"varForSampling"],
                 VCOV_sampling = VCOV_sampling,
                 VCOV_vT = VCOV_vT,
                 trace = FALSE)   
    r <- do.call("rbind", r$selinf)
    r$variable <- extractSelFun(selected_model)[which]
    r$sel = paste(selection, collapse=" + ")
    r$naive_p = summary(selected_model)$coefficients[r$variable,5]
    return(r)
    
  }, mc.cores = nrCores)
  
  
  res <- res[!sapply(res, is.null)]
  res <- do.call("rbind", res)
  if(saveResults) saveRDS(res, file=paste0(saveName,i,".RDS"))
  
  # res <- do.call("rbind", lapply(res, function(r){
  #   
  #   sel <- r$sel
  #   r <- r$selinf
  #   nr <- names(r)
  #   r <- do.call("rbind", lapply(r, function(x) x[,c("pval"),drop=F]))
  #   r$var = nr
  #   r$sel = sel
  #   rmelt = melt(r, id.vars = c("sel", "var"))[,c(1,2,4)]
  #   return(rmelt)
  #   
  # }))
  
}

### plot result in the end

if(plot_final_result){
  
  lf <- list.files(path = "sim_results/mm", pattern = saveName, full.names = T)
  lres <- lapply(lf, function(l){
    
    d <- readRDS(l)
    d$setting_nr = gsub(regexpr,"\\1",l)
    return(d)
    
  })
  
  res <- do.call("rbind", lres)
  settings$setting_nr = 1:nrow(settings)
  res <- merge(res, settings, by = "setting_nr")
  res[,c(1:6,9)] <- lapply(res[,c(1:6,9)], as.numeric)
  res[,setdiff(7:18,9)] <- lapply(res[,setdiff(7:18,9)], as.factor)
  
  print(
    gg <- ggplot(data = res,
                 aes(sample = pval,
                     colour = interaction(efficient, varForSampling),
                     linetype = bayesian)) +
      geom_abline(slope = 1, intercept = 0, linetype=2) +
      geom_qq(#,
        distribution = stats::qunif, size = 0.3, geom="line") +
      theme_bw() +# coord_flip() +
      xlab("Expected Quantiles") + ylab("Observed p-values") +
      facet_grid(sel*variable ~ SNR*marginal*sd)
  )
  
  
  saveRDS(gg, file=paste0(plotName,".RDS"))
  
}