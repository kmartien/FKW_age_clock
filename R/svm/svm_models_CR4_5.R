library(tidyverse)
library(mgcv)
library(e1071)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')
load('R/svm/svm.tuning.rda')

sites.2.use <- 'RFsites' #'All' or 'RFsites'
nrep <- 1000
ncores <- 10

svm.params <- svm.tuning$CR4_5[[sites.2.use]]$best.parameters

sites <- sites.to.keep
if(sites.2.use == 'RFsites'){
  # select important sites from Random Forest
  sites <- readRDS('R/random forest/rf_site_importance_CR4&5.rds') |> 
    filter(pval <= 0.1) |> 
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

train.df <- filter(model.df, age.confidence %in% 4:5)

# Best age and methylation estimates --------------------------------------

predictAllIDsSVM(train.df, model.df, sites, 'age.best', svm.params) |> 
  saveRDS(paste0('R/svm/svm_best_CR4_5_', sites.2.use,'.rds'))


# Random age and best methylation estimates -------------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation - only use random age
  ran.df <- model.df |> 
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |> 
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', svm.params)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/svm/svm_ran_age_CR4_5_', sites.2.use, '.rds'))


# Random age and random methylation estimates -----------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', svm.params)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/svm/svm_ran_age_meth_CR4_5_', sites.2.use, '.rds'))
