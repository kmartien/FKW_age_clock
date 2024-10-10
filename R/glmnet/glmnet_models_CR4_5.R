rm(list = ls())
library(tidyverse)
library(glmnet)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')

minCR <- 4
sites.2.use <- 'RFsites' #'Allsites', 'RFsites', 'glmnet.5', 'gamsites
age.transform <- 'ln'
weight <- 'none' # 'CR', 'inv.var', 'sn.wt', 'none'
if (weight == 'ci.wt') age.df$ci.wt <- calc.ci.wt(age.df)
nrep <- 1000
ncores <- 10

optimum.alpha <- readRDS('R/glmnet/optim.alpha.rds')$minCR4[[sites.2.use]]
sites <- sites.to.keep
if(sites.2.use != 'Allsites') sites <- selectCpGsites(sites.2.use)

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)|> 
  mutate(
    wt = if(weight == 'ci.wt') ci.wt else {
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

train.df <- filter(model.df, age.confidence >= minCR)

# Best age and methylation estimates --------------------------------------

# predictAllIDsENR(train.df, model.df, sites, 'age.best', optimum.alpha, age.transform) |>
#  saveRDS(paste0('R/glmnet/glmnet_best_minCR', minCR, '_', sites.2.use,'_', age.transform, '_', weight, '.rds'))


# Random age and best methylation estimates -------------------------------

pred.ranAge <- parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation - only use random age
  ran.df <- model.df |>
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |>
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )

  train.df <- filter(ran.df, age.confidence >= minCR)
  predictAllIDsENR(train.df, ran.df, sites, 'age.ran', optimum.alpha, age.transform)
}, mc.cores = ncores) |>
  bind_rows()
saveRDS(pred.ranAge, paste0('R/glmnet/glmnet_ranAge_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))


# Random age and random methylation estimates -----------------------------

pred.ranAgeMeth <- parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params)

  train.df <- filter(ran.df, age.confidence >= minCR)
  predictAllIDsENR(train.df, ran.df, sites, 'age.ran', optimum.alpha, age.transform)
}, mc.cores = ncores) |>
  bind_rows()
saveRDS(pred.ranAgeMeth, paste0('R/glmnet/glmnet_ranAgeMeth_minCR', minCR, '_', sites.2.use, '_', age.transform, '_', weight, '.rds'))
