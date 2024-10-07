rm(list = ls())
library(tidyverse)
library(mgcv)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")

sites.2.use <- 'RFsites' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites'
minCR <- 4
age.transform <- 'ln'
weight <- 'none' # 'CR', 'inv.var', 'sn.wt', 'none'
nrep <- 1000
ncores <- 10

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep) |> 
  mutate(
    wt = if(weight == 'inv.var') 1/age.var else {
      if (weight == 'CR') age.confidence else {
        if (weight == 'sn.wt') confidence.wt else 1
      }
    }) 

# form model data with mean logit(Pr(methylation))
model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  )

train.df <- filter(model.df, age.confidence >= minCR)

sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

# Best age and methylation estimates --------------------------------------

# message(format(Sys.time()), ' : Best - All')
# pred <- if (minCR == 2){
#   # LOO cross validation for all samples
#   lapply(model.df$swfsc.id, function(cv.id) {
#     fitTrainGAM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', age.transform) |> 
#       predictTestGAM(filter(model.df, swfsc.id == cv.id), 'age.best', age.transform)
#   })|> bind_rows() 
# } else {predictAllIDsGAM(train.df, model.df, sites, 'age.best', age.transform)}
# saveRDS(pred, paste0('R/gam/gam_best_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))

# Random age and best methylation estimates -------------------------------

message(format(Sys.time()), ' : RanAge - All')
pred <- parallel::mclapply(1:nrep, function(j) {
  ran.df <- model.df |>
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |>
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )
  if (minCR == 2){
    # LOO cross validation for all samples
    lapply(model.df$swfsc.id, function(cv.id) {
      fitTrainGAM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', age.transform) |> 
        predictTestGAM(filter(ran.df, swfsc.id == cv.id), 'age.ran', age.transform)
    })|> bind_rows() 
  } else {predictAllIDsGAM(filter(ran.df, age.confidence >= minCR), ran.df, sites, 'age.ran', age.transform)}
}, mc.cores = ncores) |> bind_rows()
saveRDS(pred, paste0('R/gam/gam_ranAge_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))

# Random age and random methylation estimates -----------------------------

message(format(Sys.time()), ' : RanAgeMeth - All')
pred <- parallel::mclapply(1:nrep, function(j) {
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)
  if (minCR == 2){
    # LOO cross validation for all samples
    lapply(model.df$swfsc.id, function(cv.id) {
      fitTrainGAM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', age.transform) |> 
        predictTestGAM(filter(ran.df, swfsc.id == cv.id), 'age.ran', age.transform)
    })|> bind_rows() 
  } else {predictAllIDsGAM(filter(ran.df, age.confidence >= minCR), ran.df, sites, 'age.ran', age.transform)}
}, mc.cores = ncores) |> bind_rows()
saveRDS(pred, paste0('R/gam/gam_ranAgeMeth_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
