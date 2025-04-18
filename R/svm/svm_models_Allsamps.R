rm(list=ls())
library(tidyverse)
library(mgcv)
library(e1071)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')
load('R/svm/svm.tuning.rda')

minCR <- 2

sites.2.use <- 'RFsites' #'Allsites' or 'RFsites', 'glmnet.5', 'gamsites'
age.transform <- 'ln'
weight <- 'none'
nrep <- 1000
ncores <- 10

svm.params <- svm.tuning$Allsamps[[sites.2.use]]$best.parameters

sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)|> 
  mutate(wt = 1) 

model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) 

# Best age and methylation estimates --------------------------------------

##########################################
# Loop through all site selection options
lapply(c('Allsites', 'RFsites', 'glmnet.5'), function(sites.2.use){
  sites <- sites.to.keep
  if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)
# LOO cross validation for all samples
pred.best <- lapply(model.df$swfsc.id, function(cv.id) {
  fitTrainSVM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', svm.params, age.transform) |>
    predictTestSVM(filter(model.df, swfsc.id == cv.id), sites, 'age.best', age.transform)
}) |> bind_rows()
saveRDS(pred.best, paste0('R/svm/svm_best_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
})

# Random age and best methylation estimates -------------------------------

pred.ranAge <- parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation - only use random age
  ran.df <- model.df |>
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |>
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )
  
  lapply(ran.df$swfsc.id, function(cv.id) {
    fitTrainSVM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', svm.params, age.transform) |>
      predictTestSVM(filter(ran.df, swfsc.id == cv.id), sites, 'age.ran', age.transform)
  }) |> bind_rows()
}, mc.cores = ncores) |>
  bind_rows()
saveRDS(pred.ranAge, paste0('R/svm/svm_ranAge_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))


# Random age and random methylation estimates -----------------------------

pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)
  
  lapply(ran.df$swfsc.id, function(cv.id) {
    fitTrainSVM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', svm.params, age.transform) |>
      predictTestSVM(filter(ran.df, swfsc.id == cv.id), sites, 'age.ran', age.transform)
  }) |> bind_rows()
}, mc.cores = ncores) |>
  bind_rows()
saveRDS(pred.ranAgeMeth, paste0('R/svm/svm_ranAgeMeth_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))

