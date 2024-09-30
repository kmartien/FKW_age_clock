library(tidyverse)
library(mgcv)
library(e1071)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

sites <- sites.to.keep

model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) 

sites.2.use <- list("Allsites", "RFsites", 'glmnet.5', 'gamsites')

svm.tuning <- lapply(c(2, 4), function(cr){
  print(paste0('minCR', cr))
  
  #select training samples
  train.df <- filter(model.df, age.confidence >= cr)
  samps <- ifelse(cr == 2, "Allsamps", "CR4&5")
  
  t <- lapply(sites.2.use, function(s){
    print(s)
    sites <- sites.to.keep
    if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)
    tune.obj <- tune(svm,
         age.best ~ .,
         data = select(train.df, c(age.best, all_of(sites))),
         ranges = list(
           cost = 10^(seq(-4, 5, 0.1)),
           gamma = 10^(seq(-5, 4, 0.1))),
         tunecontrol = tune.control(sampling = "cross"),
         cross = 10)
    return(list(best.parameters = tune.obj$best.parameters, performances = tune.obj$performances))
  })
  names(t) <- sites.2.use
  return(t)
})
names(svm.tuning) <- c("Allsamps", "CR4_5")
save(svm.tuning, file = 'R/svm/svm.tuning.rda')
