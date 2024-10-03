rm(list=ls())
library(tidyverse)
library(randomForest)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')

minCR <- 2
sites.2.use <- 'RFsites' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites'
age.transform <- 'ln'
weight <- 'none' # 'CR', 'inv.var', 'sn.wt', 'none'
nrep <- 1000
ncores <- 10

sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

load(paste0('R/rf_tuning/rf_optim_params_minCR', minCR, '.rda'))
rf.params <- rf.params[[sites.2.use]]

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)|> 
  mutate(
    wt = if(weight == 'inv.var') 1/age.var else {
      if (weight == 'CR') age.confidence else {
        if (weight == 'sn.wt') confidence.wt else 1
      }
    }) 

model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  )


# Best age and methylation estimates --------------------------------------

print('Best')
train.df <- filter(model.df, age.confidence >= minCR)
pred <- if(minCR == 2) {
  # OOB predictions for training samples
  predictTestRF(fit = NULL, train.df, sites, 'age.best', age.transform)
} else {
  predictAllIDsRF(train.df, model.df, sites, 'age.best', rf.params, age.transform)  
  }
saveRDS(pred, paste0('R/rf/rf_best_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))


# Random age and best methylation estimates -------------------------------

print('RanAge')
parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation - only use random age
  ran.df <- model.df |>
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |>
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )

  train.df <- filter(ran.df, age.confidence >= minCR)
  if(minCR == 2) {
    # OOB predictions for training samples
    predictTestRF(fit = NULL, train.df, sites, 'age.ran', age.transform)
  } else {
    predictAllIDsRF(train.df, ran.df, sites, 'age.ran', rf.params, age.transform)  
  }
}, mc.cores = ncores) |>
  bind_rows() |>
  saveRDS(paste0('R/rf/rf_ranAge_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))


# Random age and random methylation estimates -----------------------------

print('RanAgeMeth')
parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)

  train.df <- filter(ran.df, age.confidence >= minCR)
  if(minCR == 2) {
    # OOB predictions for training samples
    predictTestRF(fit = NULL, train.df, sites, 'age.ran', age.transform)
  } else {
    predictAllIDsRF(train.df, ran.df, sites, 'age.ran', rf.params, age.transform)  
  }
}, mc.cores = ncores) |>
  bind_rows() |>
  saveRDS(paste0('R/rf/rf_ranAgeMeth_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
