rm(list = ls())
library(tidyverse)
library(mgcv)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")

sites.2.use <- 'glmnet.5' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites'
minCR <- 4
age.transform <- 'ln'
weight <- 'none' # 'CR', 'inv.var', 'sn.wt', 'none'
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

train.df <- filter(model.df, age.confidence >= minCR)

sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

# best age and methylation
message(format(Sys.time()), ' : Best - All')
pred <- if (minCR == 2){
  # LOO cross validation for all samples
  lapply(model.df$swfsc.id, function(cv.id) {
    fitTrainGAM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', age.transform) |> 
      predictTestGAM(filter(model.df, swfsc.id == cv.id), 'age.best', age.transform)
  })|> bind_rows() 
} else {predictAllIDsGAM(train.df, model.df, sites, 'age.best', age.transform)}
saveRDS(pred, paste0('R/gam/gam_best_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))

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
