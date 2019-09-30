# Load libraries and source functions
# TODO: Clean up and use library for selfmade
library(lme4)
library(lmerTest)
# library(cAIC4)
# library(gamm4)
library(Matrix)
library(MASS)
library(ggplot2)
library(dplyr)
library(parallel)
library(reshape2)
library(mvtnorm)
library(selfmade)

set.seed(0) 

nrSamples = 500
nrSimIter = 100
saveResults = TRUE


n <- 150
p <- 6
nrsubj <- 5

X <- cbind(rep(1,n), matrix(rnorm(n*(p)), ncol = p))

beta <- c(1, 2, -1, -2, rep(0, p-3))
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
plotName = "ggres_check_mixed_model"
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
settings <- settings[settings$bayesian==FALSE,]
settings <- settings[settings$efficient==TRUE,]
settings <- unique(settings)


# set factors to character
settings[,sapply(settings, class)=="factor"] <- 
  lapply(settings[,sapply(settings, class)=="factor"],
         as.character)
VCOV_vT = VCOV_sampling = NULL
settings$sd = signal / settings$SNR

####################################################################
########## Generate models with false positive selection ###########
####################################################################

tau <- settings$tau[1]
covTaus <- 0.5*(0.5*tau) # = correlation * sqrt(tau * tauSlope)

vcovTRUE <- matrix(c(tau, covTaus, covTaus, 0.5 * tau), ncol = 2)

form <- as.formula(
  c("~ V1 + V2 + V3 + V4 + V5 + V6 + (1 + V2|subjind)")
)

modFun <- function(yy) 
{
  dat$y <- as.numeric(yy)
  lmer(y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + 
         (1 + V2|subjind), REML = FALSE, data = dat)
}

selFun <- function(mod)
{
  
  suppressWarnings(suppressMessages(attr(
    # use lmerTest:::step.lmerModLmerTest directly
    # as the packages overloading of step
    # does not work for the latest version
    lmerTest:::step.lmerModLmerTest(mod, reduce.random = FALSE), "model")))
  
}

extractSelFun <- function(this_mod){
  
  if(class(this_mod)=="lm") 
    return(attr(this_mod$coefficients, "names")) else
      return(c(names(fixef(this_mod)), names(getME(this_mod, "theta"))))
  
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
    return(0)
    
  }
  
}

ys_interest = vector("list", length(unique(settings$sd)))
ys_correct = vector("list", length(unique(settings$sd)))
sel_interest = vector("list", length(unique(settings$sd)))
sel_correct = vector("list", length(unique(settings$sd)))
which_interest = vector("list", length(unique(settings$sd)))
which_correct = vector("list", length(unique(settings$sd)))
mod_interest = vector("list", length(unique(settings$sd)))
mod_correct = vector("list", length(unique(settings$sd)))

k = 1
totalmult = 1000

for(sigma in unique(settings$sd)){
  
  print(paste0("sigma: ", round(sigma,2)))
  total <- nrSimIter
  m = 1
  
  cholVCT <- chol(vcovTRUE)/sigma
  theta <- c(cholVCT)[-2]
  names(theta)<-c("subjind.(Intercept)",
                  "subjind.V2.(Intercept)",
                  "subjind.V2")
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  set.seed(m)
  this_ys <- suppressWarnings(suppressMessages(
    simulate(form, newdata = dat, family = "gaussian",
             newparams = list(theta = theta, beta = beta, 
                              sigma = sigma),
             nsim = totalmult*total)
  ))
  ys_interest[[k]] <- which_interest[[k]] <- 
    ys_correct[[k]] <- which_correct[[k]] <- 
    sel_interest[[k]] <- sel_correct[[k]] <- 
    mod_correct[[k]] <- mod_interest[[k]] <- list()
  
  j <- 1
  while(m < total + 1 & j < totalmult*total + 1){
    
    # check selected model
    this_y = this_ys[,j]
    mod = selFun(modFun(this_y))
    if(length(mod@optinfo$conv$lme4)>0){
      j = j + 1
      next
    }
    selection = extractSelFun(mod)  
    which = suppressWarnings(nameFun(mod))
    if(which==5){
      
      ys_interest[[k]] = c(ys_interest[[k]], list(j))
      which_interest[[k]] = c(which_interest[[k]], list(which))
      sel_interest[[k]] = c(sel_interest[[k]], list(selection))
      mod_interest[[k]] = c(mod_interest[[k]], mod)
      setTxtProgressBar(pb, m)
      m = m + 1
      
    }else if(which==2 & length(ys_correct[[k]]) < 101){
      
      ys_correct[[k]] = c(ys_correct[[k]], list(j))
      which_correct[[k]] = c(which_correct[[k]], list(which))
      sel_correct[[k]] = c(sel_correct[[k]], list(selection))
      mod_correct[[k]] = c(mod_correct[[k]], mod)
      
    }
    j = j + 1


  }
  close(pb)
  ys_interest[[k]] = this_ys[,c(unlist(ys_interest[[k]]), 
                                unlist(ys_correct[[k]])[1:50])]
  which_interest[[k]] = c(unlist(which_interest[[k]]), 
                          unlist(which_correct[[k]])[1:50])
  sel_interest[[k]] = c(sel_interest[[k]], sel_correct[[k]][1:50])
  mod_interest[[k]] = c(mod_interest[[k]], mod_correct[[k]][1:50])
  k = k + 1
  
}
    
####################################################################
############### Start iterating over simulations ###################
####################################################################

settings$y_mapping <- sapply(settings$sd, function(x) 
  which(unique(settings$sd)==x))

for(i in 1:nrow(settings)){
  
  if(settings[i,"varForSampling"]=="supplied"){
    # supply true value
    VCOV_sampling = settings[i,"sd"]^2 * diag(rep(1,n))
    
  }
  
  ys_list = ys_interest[[settings$y_mapping[i]]]
  which_list = which_interest[[settings$y_mapping[i]]]
  # added the following for testing different true effects
  which_list[which_list==2] <- rep(2:3, each=sum(which_list==2)/2)
  sel_list = sel_interest[[settings$y_mapping[i]]]
  mod_list = mod_interest[[settings$y_mapping[i]]]
  
  
  res <- mclapply(1:length(ys_list), function(j){
    
    st <- Sys.time()
    sigma <- settings$sd[i]
    
    ### define selection functions
    modFun <- function(yy) 
    {
      dat$y <- as.numeric(yy)
      lmer(y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + 
             (1 + V2|subjind), REML = FALSE, data = dat)
    }
    
    selFun <- function(mod)
    {
      
      attr(lmerTest:::step.lmerModLmerTest(mod, reduce.random = FALSE), "model")
      
    }
    
    extractSelFun <- function(this_mod){
      
      if(class(this_mod)=="lm") 
        return(attr(this_mod$coefficients, "names")) else
          return(c(names(fixef(this_mod)), names(getME(this_mod, "theta"))))
      
    }
    
    selection = sel_list[[j]]
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
        if(length(sel_vars)==4) return(1) else
          return(5)
        
      }else{
        
        # not correct model selected
        if(length(sel_vars)>1) return(2) else return(1)
        
      }
      
    }

    r <- mocasin(mod = mod_list[[j]], 
                 this_y = ys_list[[j]],
                 checkFun = checkFun, 
                 nrSamples = nrSamples,
                 which = which_list[[j]],
                 bayesian = settings[i,"bayesian"],
                 conditional = settings[i,"conditional"],
                 efficient = settings[i,"efficient"],
                 varForSampling = settings[i,"varForSampling"],
                 VCOV_sampling = VCOV_sampling,
                 VCOV_vT = VCOV_vT,
                 trace = FALSE)   
    r <- do.call("rbind", r$selinf)
    r$variable <- sel_list[[j]][which_list[[j]]]
    r$sel = paste(sel_list[[j]], collapse=" + ")
    r$naive_p = summary(mod_list[[j]])$coefficients[r$variable,5]
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
  
  find_true_setting <- function(x)
  {
    
    fixef_pres = all(sapply(paste0("V",1:3), function(var) grepl(pattern = var, 
                                                                 x, fixed = TRUE)))
    ranef_pres = grepl("subjind.(Intercept)", x, fixed = T) & 
      grepl("subjind.V2.(Intercept)", x, fixed = T) & 
      grepl("subjind.V2", x, fixed = T)
    
    return(c(fixef_pres, ranef_pres))
    
  }
  
  seltype = as.data.frame(t(sapply(res$sel, find_true_setting)))
  colnames(seltype) = c("allfixed", "allrandom")
  res = cbind(res, seltype)
  
  res <- res[res$allfixed & res$allrandom,]
  res$variable <- as.character(droplevels(res$variable))
  res$variable[res$variable%in%paste0("V",4:6)] <- "noise"
  res$variable <- as.factor(res$variable)
  
  # print(
  #   gg1 <- ggplot(data = res[res$bayesian=="TRUE",],
  #                 aes(sample = pval,
  #                     colour = varForSampling
  #                     #linetype = bayesian
  #                 )) +
  #     geom_abline(slope = 1, intercept = 0, linetype=2) +
  #     geom_qq(#,
  #       distribution = stats::qunif, size = 0.8, geom="line") +
  #     theme_bw() +# coord_flip() +
  #     xlab("Expected Quantiles") + ylab("Observed p-values") +
  #     facet_grid(variable ~ varForSampling*SNR)
  # )
  
  levels(res$variable) <- c("noise variables",
                            "signal variable 1", 
                            "signal variable 2")
  res$SNRfac <- factor(res$SNR, levels = unique(res$SNR),
                       labels = paste0("SNR = ", unique(res$SNR)))
  levels(res$marginal) <- c("conditional", "marginal")
  levels(res$varForSampling) <- c("Model Estimate",
                                  "ICM",
                                  "Truth",
                                  "Var(Y)")
  
  mres = melt(res %>% 
                select(pval, variable, naive_p,
                       varForSampling, marginal,
                       SNRfac), 
              id.vars=c("variable", "varForSampling",
                        "marginal", "SNRfac"), variable.name = "p_val_type")
  
  levels(mres$p_val_type) <- c("selective", "naive")
  
  print(
    gg2 <- ggplot(data = mres,
                  aes(sample = value,
                      colour = varForSampling,
                      linetype = p_val_type
                  )) +
      geom_abline(slope = 1, intercept = 0, linetype=2) +
      geom_qq(#,
        distribution = stats::qunif, size = 0.8, geom="line") +
      theme_bw() +# coord_flip() +
      xlab("Expected Quantiles") + ylab("Observed p-values") +
      facet_grid(variable ~ SNRfac*marginal) + 
      labs(colour = "Variance", linetype="p-value type")
  )
  
  saveRDS(gg2, file=paste0(plotName,".RDS"))
  
}
