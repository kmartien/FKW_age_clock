rm(list = ls())
library(tidyverse)
library(mgcv)
source('../misc_funcs.R')
load("../../data/age_and_methylation_data.rdata")

# select important sites from Random Forest
sites <- readRDS('../random forest/rf_site_importance.rds') |> 
  filter(pval <= 0.1) |> 
  pull('loc.site')

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

# form model data with mean logit(Pr(methylation))
model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) 

nrep <- 100
ncores <- 10

# Best age and methylation estimates --------------------------------------

train.df <- filter(model.df, age.confidence %in% 4:5)
predictAllIDsGAM(train.df, model.df, sites, 'age.best') |> 
  saveRDS('gam_best.rds')


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
  predictAllIDsGAM(train.df, ran.df, sites, 'age.ran')
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS('gam_ran_age.rds')


# Random age and random methylation estimates -----------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsGAM(train.df, ran.df, sites, 'age.ran')
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS('gam_ran_age_meth.rds')
