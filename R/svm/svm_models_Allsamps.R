library(tidyverse)
library(mgcv)
library(e1071)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')
load('R/svm/svm.tuning.rda')

sites.2.use <- 'RFsites' #'Allsites' or 'RFsites'
nrep <- 1000
ncores <- 10

svm.params <- svm.tuning$Allsamps[[sites.2.use]]$best.parameters

sites <- sites.to.keep
if(sites.2.use == 'RFsites'){
  # select important sites from Random Forest
  sites <- readRDS('R/rf_tuning/rf_site_importance_Allsamps.rds') |> 
    filter(pval <= 0.05) |> 
    pull('loc.site')
}

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)
model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) 

# Best age and methylation estimates --------------------------------------

# LOO cross validation for all samples
lapply(model.df$swfsc.id, function(cv.id) {
  fitTrainSVM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', svm.params) |> 
    predictTestSVM(filter(model.df, swfsc.id == cv.id), sites, 'age.best')
}) |> bind_rows() |>
  saveRDS(paste0('R/svm/svm_best_Allsamps_', sites.2.use,'.rds'))


# Random age and best methylation estimates -------------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation - only use random age
  ran.df <- model.df |> 
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |> 
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )
  
  lapply(ran.df$swfsc.id, function(cv.id) {
    fitTrainSVM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', svm.params) |> 
      predictTestSVM(filter(ran.df, swfsc.id == cv.id), sites, 'age.ran')
  }) |> bind_rows()
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/svm/svm_ran_age_Allsamps_', sites.2.use, '.rds'))


# Random age and random methylation estimates -----------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
  
  lapply(ran.df$swfsc.id, function(cv.id) {
    fitTrainSVM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', svm.params) |> 
      predictTestSVM(filter(ran.df, swfsc.id == cv.id), sites, 'age.ran')
  }) |> bind_rows()
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/svm/svm_ran_age_meth_Allsamps_', sites.2.use, '.rds'))
