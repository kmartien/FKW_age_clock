library(tidyverse)
library(mgcv)
library(e1071)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')
load('R/svm/svm.tuning.rda')

sites.2.use <- 'glmnet.5' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites'
age.transform <- 'ln'
weight <- 'none'
nrep <- 1000
ncores <- 10

# NEED TO FIX THE NEXT TO LINE FOR SVM TUNING
svm.params <- svm.tuning$CR4_5$RFsites$best.parameters

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
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  )


train.df <- filter(model.df, age.confidence %in% 4:5)

# Best age and methylation estimates --------------------------------------

predictAllIDsSVM(train.df, model.df, sites, 'age.best', svm.params, age.transform) |> 
  saveRDS(paste0('R/svm/svm_best_minCR4_', sites.2.use,'_', age.transform, '.rds'))


# # Random age and best methylation estimates -------------------------------
# 
# parallel::mclapply(1:nrep, function(j) {
#   # random sample of ages and methylation - only use random age
#   ran.df <- model.df |> 
#     left_join(
#       sampleAgeMeth(age.df, logit.meth.normal.params) |> 
#         select(swfsc.id, age.ran),
#       by = 'swfsc.id'
#     )
#   
#   train.df <- filter(ran.df, age.confidence %in% 4:5)
#   predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', svm.params)
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS(paste0('R/svm/svm_ran_age_minCR4_', sites.2.use, '.rds'))
# 
# 
# # Random age and random methylation estimates -----------------------------
# 
# parallel::mclapply(1:nrep, function(j) {
#   # random sample of ages and methylation
#   ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
#   
#   train.df <- filter(ran.df, age.confidence %in% 4:5)
#   predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', svm.params)
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS(paste0('R/svm/svm_ran_age_meth_minCR4_', sites.2.use, '.rds'))
