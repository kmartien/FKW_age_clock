rm(list = ls())
library(tidyverse)
library(mgcv)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")

sites.2.use <- 'glmnet.5' #'Allsites', 'RFsites', 'glmnet.5'
minCR <- 4
age.transform <- 'ln'
weight <- 'sn.wt' # 'CR', 'inv.var', 'sn.wt', 'none'
nrep <- 1000
ncores <- 10

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

# form model data with mean logit(Pr(methylation))
model.df <- age.df |> 
  mutate(
    wt = if(weight == 'inv.var') 1/age.var else {
      if (weight == 'CR') age.confidence else {
        if (weight == 'sn.wt') confidence.wt else 1
      }
    }) |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  )

train.df <- filter(model.df, age.confidence %in% 4:5)

# Using sites from RF trained on all samples ------------------------------

sites <- sites.to.keep
if(sites.2.use == 'RFsites'){
  # select important sites from Random Forest tuned with all samples
  sites <- readRDS('R/rf_tuning/rf_site_importance_Allsamps.rds') |>
    filter(pval <= 0.05) |>
    pull('loc.site')
}
if(sites.2.use == 'glmnet.5'){
  # select chosen sites from glmnet tuned with all samples, alpha = 0.5
  sites <- readRDS('R/glmnet/glmnet.chosen.sites.rds')$alpha.5$minCR2
}
if(sites.2.use == 'gamsites'){
  # select chosen sites from gam.by.site
  sites <- readRDS('R/gam/gam_significant_sites.rds') |>
    filter(r.sq >= 0.35) |>
    pull('loc.site')
}

# best age and methylation
message(format(Sys.time()), ' : Best - All')
predictAllIDsGAM(train.df, model.df, sites, 'age.best', age.transform) |>
  saveRDS(paste0('R/gam/gam_best_minCR4_', sites.2.use, '_', age.transform, '_', weight, '.rds'))

# # random age and best methylation
# message(format(Sys.time()), ' : RanAge - All')
# parallel::mclapply(1:nrep, function(j) {
#   ran.df <- model.df |> 
#     left_join(
#       sampleAgeMeth(age.df, logit.meth.normal.params) |> 
#         select(swfsc.id, age.ran),
#       by = 'swfsc.id'
#     ) 
#   predictAllIDsGAM(filter(ran.df, age.confidence %in% 4:5), ran.df, sites_Allsamps, 'age.ran')
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS('gam_ran_age_Allsamps.rds')
# 
# # random age and random methylation
# message(format(Sys.time()), ' : RanAgeMeth - All')
# parallel::mclapply(1:nrep, function(j) {
#   ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
#   predictAllIDsGAM(filter(ran.df, age.confidence %in% 4:5), ran.df, sites_Allsamps, 'age.ran')
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS('gam_ran_age_meth_Allsamps.rds')
# 
# 
# 
# 
# # Using sites from RF trained on CR 4 & 5 samples -------------------------
# 
# # select important sites from Random Forest trained on CR 4 & 5 samples
# sites_CR4and5 <- readRDS('../random forest/rf_site_importance_CR4&5.rds') |> 
#   filter(pval <= 0.05) |> 
#   pull('loc.site')
# 
# # best age and methylation
# message(format(Sys.time()), ' : Best - CR4and5')
# predictAllIDsGAM(train.df, model.df, sites_CR4and5, 'age.best') |> 
#   saveRDS('gam_best_CR4and5.rds')
# 
# # random age and best methylation
# message(format(Sys.time()), ' : RanAge - CR4and5')
# parallel::mclapply(1:nrep, function(j) {
#   ran.df <- model.df |> 
#     left_join(
#       sampleAgeMeth(age.df, logit.meth.normal.params) |> 
#         select(swfsc.id, age.ran),
#       by = 'swfsc.id'
#     )
#   predictAllIDsGAM(filter(ran.df, age.confidence %in% 4:5), ran.df, sites_CR4and5, 'age.ran')
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS('gam_ran_age_CR4and5.rds')
# 
# # random age and random methylation
# message(format(Sys.time()), ' : RanAgeMeth - CR4and5')
# parallel::mclapply(1:nrep, function(j) {
#   ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
#   predictAllIDsGAM(filter(ran.df, age.confidence %in% 4:5), ran.df, sites_CR4and5, 'age.ran')
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS('gam_ran_age_meth_CR4and5.rds')
