rm(list = ls())
library(tidyverse)
library(mgcv)
library(randomForest)
source('R/misc_funcs.R')
load("R/data/age_and_methylation_data.rdata")

minCR <- 4

sites.2.use <- "All" #"All" or "RFsites"
load(paste0('R/random forest/rf_optim_params_CR4&5_', sites.2.use, '.rda'))
#rf.params <- tibble(mtry = floor(length(sites)/3), sampsize = floor(50 * 0.632)) 

sites <- sites.to.keep
if(sites.2.use == "RFsites"){
  # select important sites from Random Forest
  sites <- readRDS('R/random forest/rf_site_importance.rds') |> 
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

nrep <- 1000
ncores <- 10

# Best age and methylation estimates --------------------------------------

train.df <- filter(model.df, age.confidence <= minCR)
predictAllIDsRF(train.df, model.df, sites, 'age.best', rf.params) |> 
  saveRDS(paste0('R/random forest/rf_best_minCR', minCR, '_', sites.2.use, '.rds'))


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
  predictAllIDsRF(train.df, ran.df, sites, 'age.ran', rf.params)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/random forest/rf_ran_age_minCR', minCR, '_', sites.2.use, '.rds'))


# Random age and random methylation estimates -----------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsRF(train.df, ran.df, sites, 'age.ran', rf.params)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/random forest/rf_ran_age_meth_minCR', minCR, '_', sites.2.use, '.rds'))
