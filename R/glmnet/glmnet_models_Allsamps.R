rm(list=ls())
library(tidyverse)
library(mgcv)
library(glmnet)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')

minCR <- 2
sites.2.use <- 'RFsites' #'Allsites' or 'RFsites'
age.transform <- 'ln'
nrep <- 1000
ncores <- 10

optimum.alpha <- readRDS('R/glmnet/optimum.alpha.rds')$minCR2
sites <- sites.to.keep
if(sites.2.use == 'RFsites'){
  # select important sites from Random Forest
  sites <- readRDS('R/glmnet/glmnet.chosen.sites.rds')$minCR2
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

lapply(model.df$swfsc.id, function(cv.id) {
  fitTrainENR(filter(model.df, swfsc.id != cv.id), sites, 'age.best', optimum.alpha$par, age.transform) |> 
    predictTestENR(filter(model.df, swfsc.id == cv.id), sites, 'age.best', age.transform)
}) |> 
  bind_rows() |>
  saveRDS(paste0('R/glmnet/glmnet_best_minCR', minCR, '_', sites.2.use,'_', age.transform, '.rds'))


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
#   train.df <- filter(ran.df, age.confidence >= minCR)
#   lapply(ran.df$swfsc.id, function(cv.id) {
#     fitTrainENR(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', optimum.alpha$par) |> 
#       predictTestENR(filter(ran.df, swfsc.id == cv.id), sites, 'age.ran')
#   }) |> 
#     bind_rows() 
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS(paste0('R/glmnet/glmnet_ran_age_minCR', minCR, '_', sites.2.use, '.rds'))
# 
# 
# # Random age and random methylation estimates -----------------------------
# 
# parallel::mclapply(1:nrep, function(j) {
#   # random sample of ages and methylation
#   ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
#   
#   train.df <- filter(ran.df, age.confidence >= minCR)
#   lapply(ran.df$swfsc.id, function(cv.id) {
#     fitTrainENR(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', optimum.alpha$par) |> 
#       predictTestENR(filter(ran.df, swfsc.id == cv.id), sites, 'age.ran')
#   }) |> 
#     bind_rows() 
# }, mc.cores = ncores) |> 
#   bind_rows() |> 
#   saveRDS(paste0('R/glmnet/glmnet_ran_age_meth_minCR', minCR, '_', sites.2.use, '.rds'))
