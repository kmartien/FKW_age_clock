rm(list = ls())
library(tidyverse)
library(mgcv)
source('R/misc_funcs.R')
source('R/calc.ci.wt.R')
load("data/age_and_methylation_data.rdata")

sites.2.use <- 'RFsites' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites', 'glmnet.June'
minCR <- 2
age.transform <- 'ln'
weight <- 'CR' # 'CR', 'ci.wt', 'none'
if (weight == 'ci.wt') age.df$ci.wt <- calc.ci.wt(age.df)
nrep <- 1000
ncores <- 10

sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep) |> 
  mutate(
    wt = if(weight == 'ci.wt') ci.wt else {
      if (weight == 'CR') age.confidence else {
        if (weight == 'sn.wt') confidence.wt else 1
      }
    }) 

# form model data with mean logit(Pr(methylation))
model.df <- age.df |> 
  left_join(
    #  logit.meth |> 
    #    select(swfsc.id, loc.site, logit.meth) |>
    #    pivot_wider(names_from = 'loc.site', values_from = 'logit.meth'),
    # by = 'swfsc.id'
    rownames_to_column(logit.meth, var = 'swfsc.id') |> 
      select(c(swfsc.id, all_of(sites)))
    
  )
# NAs.by.sample <- sapply(1:nrow(model.df), function(s){
#   length(which(is.na(model.df[s,all_of(sites.to.keep)])))
# })
# model.df <- model.df[-which(NAs.by.sample>0),]

train.df <- filter(model.df, age.confidence >= minCR)

# Best age and methylation estimates --------------------------------------

# message(format(Sys.time()), ' : Best - All')
pred.best <- if (minCR == 2){
  # LOO cross validation for all samples
  lapply(model.df$swfsc.id, function(cv.id) {
    fitTrainGAM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', age.transform) |>
      predictTestGAM(filter(model.df, swfsc.id == cv.id), 'age.best', age.transform)
  })|> bind_rows()
} else {predictAllIDsGAM(train.df, model.df, sites, 'age.best', age.transform)}
saveRDS(pred.best, paste0('R/gam/gam_best_minCR', minCR, '_', sites.2.use, '_', age.transform, 'KKMmeth_', weight, '.rds'))

# Random age and best methylation estimates -------------------------------

# message(format(Sys.time()), ' : RanAge - All')
# pred.ranAge <- parallel::mclapply(1:nrep, function(j) {
#   ran.df <- model.df |>
#     left_join(
#       sampleAgeMeth(age.df, logit.meth.normal.params) |>
#         select(swfsc.id, age.ran),
#       by = 'swfsc.id'
#     )
#   if (minCR == 2){
#     # LOO cross validation for all samples
#     lapply(model.df$swfsc.id, function(cv.id) {
#       fitTrainGAM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', age.transform) |>
#         predictTestGAM(filter(ran.df, swfsc.id == cv.id), 'age.ran', age.transform)
#     })|> bind_rows()
#   } else {predictAllIDsGAM(filter(ran.df, age.confidence >= minCR), ran.df, sites, 'age.ran', age.transform)}
# }, mc.cores = ncores) |> bind_rows()
# saveRDS(pred.ranAge, paste0('R/gam/gam_ranAge_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
# 
# # # Random age and random methylation estimates -----------------------------
# # 
# message(format(Sys.time()), ' : RanAgeMeth - All')
# pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
#   print(j)
#   ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)
#   if (minCR == 2){
#     # LOO cross validation for all samples
#     lapply(model.df$swfsc.id, function(cv.id) {
#       fitTrainGAM(filter(ran.df, swfsc.id != cv.id), sites, 'age.ran', age.transform) |>
#         predictTestGAM(filter(ran.df, swfsc.id == cv.id), 'age.ran', age.transform)
#     })|> bind_rows()
#   } else {predictAllIDsGAM(filter(ran.df, age.confidence >= minCR), ran.df, sites, 'age.ran', age.transform)}
# }, mc.cores = ncores) |> bind_rows()
# saveRDS(pred.ranAgeMeth, paste0('R/gam/gam_ranAgeMeth_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
