rm(list = ls())
library(tidyverse)
library(e1071)
source('R/misc_funcs.R')
load("data/age_and_methylation_data.rdata")
load('R/svm/svm.tuning.rda')

r.sq_threshold <- seq(0.05,0.5, 0.025)
minCR <- 4
age.transform <- 'ln'
weight <- 'none' # 'CR', 'inv.var', 'sn.wt', 'none'
nrep <- 1000
ncores <- 10

svm.params <- svm.tuning$CR4_5$RFsites$best.parameters

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

res <- lapply(r.sq_threshold, function(r){
sites <- readRDS('R/gam/gam_significant_sites.rds') |>
  filter(r.sq >= r) |>
  pull('loc.site')

# best age and methylation
message(format(Sys.time()), ' : Best - All')
pred <- if (minCR == 2){
  # LOO cross validation for all samples
  lapply(model.df$swfsc.id, function(cv.id) {
    fitTrainSVM(filter(model.df, swfsc.id != cv.id), sites, 'age.best', svm.params, age.transform) |> 
      predictTestSVM(filter(model.df, swfsc.id == cv.id), sites, 'age.best', age.transform)
  })|> bind_rows() 
} else {predictAllIDsSVM(train.df, model.df, sites, 'age.best', svm.params, age.transform)}
pred |> mutate(threshold = r)
}) |>  bind_rows()

saveRDS(res, file = 'R/svm/test.gamsite.selection.criterion.rds')
