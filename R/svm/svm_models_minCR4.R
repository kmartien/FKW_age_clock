rm(list=ls())
library(tidyverse)
library(mgcv)
library(e1071)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')
load('data/corrected.meth.rda')
load('R/svm/svm.tuning.rda')

sites.2.use <- 'glmnet.June' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites', 'glmnet.June
minCR <- 4
age.transform <- 'ln'
weight <- 'none'
nrep <- 1000
ncores <- 10

#svm.params <- svm.tuning$CR4_5[[sites.2.use]]$best.parameters
svm.params <- svm.tuning$CR4_5$glmnet.5$best.parameters

sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

model.df <- age.df |> 
  mutate(
    wt = if(weight == 'inv.var') 1/age.var else {
      if (weight == 'CR') age.confidence else 1
  }) |> 
  left_join(
    # logit.meth.normal.params |> 
    #   select(swfsc.id, loc.site, mean.logit) |>
    #   pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    corrected.meth,
    by = 'swfsc.id'
  )
NAs.by.sample <- sapply(1:nrow(model.df), function(s){
  length(which(is.na(model.df[s,all_of(sites.to.keep)])))
})
model.df <- model.df[-which(NAs.by.sample>0),]

train.df <- filter(model.df, age.confidence >= minCR)

 # Best age and methylation estimates --------------------------------------
 
##########################################
# Loop through all site selection options
# lapply(c('glmnet.June'), function(sites.2.use){
#   sites <- sites.to.keep
#   if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)
pred <- predictAllIDsSVM(train.df, model.df, sites, 'age.best', svm.params, age.transform)
  saveRDS(pred, paste0('R/svm/svm_best_minCR', minCR, '_', sites.2.use,'_', age.transform, '_', weight, '.rds'))
#})

# Random age and best methylation estimates -------------------------------

# pred.ranAge <- parallel::mclapply(1:nrep, function(j) {
#   # random sample of ages and methylation - only use random age
#   ran.df <- model.df |>
#     left_join(
#       sampleAgeMeth(age.df, logit.meth.normal.params) |>
#         select(swfsc.id, age.ran),
#       by = 'swfsc.id'
#     )
#   
#   train.df <- filter(ran.df, age.confidence >= minCR)
#   predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', svm.params)
# }, mc.cores = ncores) |>
#   bind_rows()
# saveRDS(pred.ranAge, paste0('R/svm/svm_ranAge_minCR', minCR, '_', sites.2.use,'_', age.transform, '_', weight, '.rds'))
# 
# 
# Random age and random methylation estimates -----------------------------

 pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)

  train.df <- filter(ran.df, age.confidence >= minCR)
  predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', svm.params)
}, mc.cores = ncores) |>
  bind_rows()
saveRDS(pred.ranAgeMeth, paste0('R/svm/svm_ranAgeMeth_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
