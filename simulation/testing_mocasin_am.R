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
#library(selfmade)
library(devtools)
load_all("~/Dropbox/packages/selfmade/selfmade/R")

set.seed(0) 

n = 500
nrSamples = 500
nrSimIter = 500
saveResults = TRUE

### create data
centeredRN <- function() as.numeric(scale(rnorm(n/2), scale=F))

x1 <- c(centeredRN()%*%t(c(-1,1)))
x2 <- c(centeredRN()%*%t(c(-1,1)))
x3 <- c(centeredRN()%*%t(c(-1,1)))
x4 <- c(centeredRN()%*%t(c(-1,1)))
x5 <- c(centeredRN()%*%t(c(-1,1)))

# x6 <- rnorm(n)
dat <- data.frame(x1 = x1, x2 = x2, x3 = x3,
                  x4 = x4, x5 = x5) # , x6 = x6)
### define dgp
dgp <- function(n, sd) 1 - tanh(dat$x1) + sin(3*dat$x2) + rnorm(n=n, sd=sd)

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
plotName = "ggres_check_additive_model"
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
  varForSampling = c("est", "minMod", "varY", "supplied"),
  #varForSampling = c("est", "minMod", "varY", "supplied"),
  sd = c(1,sqrt(10))
)

settings$varInTestvec <- "est"
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
################ Finding location to be tested # ###################
####################################################################

# # generate y with almost no error
# yhere = dgp(n = n, sd = 0.05)
# dat$yhere = yhere
# # fit one model, where H_0 is true
# br3 <- gam(yhere ~ s(x1) + s(x2) + x3 + x4, data = dat)
# # predict one term
# pr = predict(br3, terms = "s(x1)", type="terms")
# # due to the fitting routine, predictions for all
# # x1 must sum to zero
# sum(pr)
# # which means that f_1(0) cannot be zero
# # instead we therefore test the value z at which 
# # sum((x1)^2 - z) == 0
# # depending on the level of noise, this z is approx. 1
# plot(pr ~ x1, pch="-")
# abline(h = 0, lty= 2)
# abline(v = 1)
# abline(v = -1)
# points(seq(-3,3,l=20),seq(-3,3,l=20)^2-1,type="l",col="red")

####################################################################
############### Start iterating over simulations ###################
####################################################################


for(i in 1:nrow(settings)){
  
  if(settings[i,"varInTestvec"]=="supplied"){
    # supply true value
    VCOV_vT = settings[i,"sd"]^2
    
  }else{
    
    VCOV_vT = NULL
    
  }
  
  if(settings[i,"varForSampling"]=="supplied"){
    # supply true value
    VCOV_sampling = settings[i,"sd"]^2 * diag(rep(1,n))
    
  }else{
    
    VCOV_sampling = NULL
    
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
                 nrlocs = c(-1,0), 
                 nrSamples = nrSamples,
                 bayesian = settings[i,"bayesian"],
                 conditional = settings[i,"conditional"],
                 efficient = settings[i,"efficient"],
                 varInTestvec = settings[i,"varInTestvec"],
                 varForSampling = settings[i,"varForSampling"],
                 VCOV_sampling = VCOV_sampling,
                 VCOV_vT = VCOV_vT,
                 trace = FALSE)   
    print(Sys.time()-st)
    r$sel = wm
    return(r)
    
  }, mc.cores = nrCores)
  
  res <- res[!sapply(res, is.null)]
  res <- do.call("rbind", lapply(res, function(r){
    
    sel <- r$sel
    nr <- names(r$vT)
    r <- r$selinf
    r <- do.call("rbind", lapply(r, function(x) x[,c("pval"),drop=F]))
    r$var = nr
    r$sel = sel
    rmelt = melt(r, id.vars = c("sel", "var"))[,c(1,2,4)]
    return(rmelt)
    
  }))
  
  if(saveResults) saveRDS(res, file=paste0(saveName,i,".RDS"))
  
}

### plot result in the end

if(plot_final_result){
  
  lf <- list.files(path = "sim_results/am/", pattern = saveName, full.names = T)
  lres <- lapply(lf, function(l){
    
    d <- readRDS(l)
    # d <- d[!sapply(d,is.null)]
    # sels <- unlist(sapply(d,"[[","sel"))
    # d <- do.call("rbind", lapply(1:length(sels), function(i){
    #   var <- names(d[[i]]$selinf)
    #   data <- do.call("rbind", d[[i]]$selinf)
    #   data$sel <- sels[i]
    #   data$var <- var
    #   return(data)
    # }))
    d$setting_nr = gsub(regexpr,"\\1",l)
    return(d)
    
  })
  
  res <- do.call("rbind", lres)
  
  settings$setting_nr = 1:nrow(settings)
  res <- merge(res, settings, by="setting_nr")
  
  res$sel <- as.factor(res$sel)
  res$varForSampling <- as.factor(res$varForSampling)
  res$varInTestvec <- as.factor(res$varInTestvec)
  res$sd = round(res$sd,2)
  
  res$sdFac = factor(res$sd, levels = unique(res$sd),
                     labels = c(expression(paste("sd = ", 1)),
                                expression(paste("sd = ", sqrt(10)))))
  res <- res[res$sel=="3",]
  res$var[res$var=="x3"] <- res$var[res$var=="x4"] <- "noise"
  res$varFac = factor(res$var, levels = unique(res$var),
                      labels = c("noise",
                                 expression(paste(x[1], "= -1")),
                                 expression(paste(x[1], "= 0")),
                                 expression(paste(x[2], "= -1")),
                                 expression(paste(x[2], "= 0"))))
  
  res$bayesian <- factor(res$bayesian, levels = c(TRUE, FALSE),
                         labels = c("Bayesian", "Classical"))
  levels(res$varInTestvec) <- c("Estimate",
                                "Truth",
                                "Var(Y)")
  levels(res$varForSampling) <- c("Estimate",
                                  "Truth",
                                  "Var(Y)")
  
  res$variances = interaction(res$varInTestvec, res$varForSampling, sep = " / ")
  
  print(
    gg <- ggplot(data = res[
      # res$varInTestvec!="Var(Y)" & 
      #   res$var%in%c("noise",
      #                "x11", "x12"),
      ],
      aes(sample = value,
                     colour = varInTestvec,
                     linetype = bayesian)) +
      geom_abline(slope = 1, intercept = 0, linetype=2) +
      geom_qq(#,
        distribution = stats::qunif, size = 0.8, geom="line") +
      theme_bw() +# coord_flip() +
      xlab("Expected Quantiles") + ylab("Observed p-values") +
      facet_grid(varFac ~ sdFac*varForSampling, labeller = label_parsed) + 
      labs(colour = "Variance used for sampling",
           colour = "Variance used for test vector") + 
      scale_linetype_manual(values = c("dashed","dotted"))
  )
  
  
  #saveRDS(gg, file=paste0(plotName,".RDS"))
  
  print(
    gg <- ggplot(data = res[
      #res$varInTestvec=="Var(Y)",
      ],
                 aes(sample = value,
                     colour = bayesian,
                     linetype = sdFac)) +
      geom_abline(slope = 1, intercept = 0, linetype=2) +
      geom_qq(#,
        distribution = stats::qunif, size = 0.8, geom="line") +
      theme_bw() +# coord_flip() +
      xlab("Expected Quantiles") + ylab("Observed p-values") +
      facet_grid(varFac ~ varForSampling, labeller = label_parsed) + 
      labs(colour = NULL) + 
      scale_linetype_manual(values = c("dashed", "dotted"),
                            labels = c(expression(paste("sd = ",1)),
                                       expression(paste("sd = ",sqrt(10))))) + 
      labs(linetype = NULL)
  )
  
  
  saveRDS(gg, file=paste0(plotName,"_2.RDS"))
  
  saveRDS(res, file="paper/results/res_check_additive_model.RDS") 
  
}