rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)
load("data/age_and_methylation_data.rdata")

minCR <- 3

#' run randomForest over sampsize and mtry grid and report deviance stats
rf.param.grid.search <- function(model.df, sites){
  # grid of sampsize and mtry to optimize MSE over
  params <- expand.grid(
    sampsize = seq(round(nrow(model.df)*.3), round(nrow(model.df)*.7), by = 1), 
    mtry = seq(round(length(sites)*.2), round(length(sites)*.5), by = 1), 
    KEEP.OUT.ATTRS = FALSE
  )
  
  parallel::mclapply(1:nrow(params), function(i){
    rf <- randomForest(
      y = model.df$age.best,
      x = model.df[, sites.to.keep],
      mtry = params$mtry[i],
      ntree = 10000,
      sampsize = params$sampsize[i],
      replace = FALSE
    )
    data.frame( 
      sampsize = params$sampsize[i],
      mtry = params$mtry[i],
      mse = rf$mse[length(rf$mse)],
      rsq = rf$rsq[length(rf$rsq)],
      pct.var = 100 * rf$rsq[length(rf$rsq)]
    )
  }, mc.cores = 6) |> 
    bind_rows() |> 
    as.data.frame()
}

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep & age.confidence >= minCR) |>
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') |> 
  mutate(inv.var = 1 / age.var)

# All sites -----------------------------------------------------------------
param.df  <- rf.param.grid.search(model.df, sites.to.keep)
save(param.df, file = paste0('R/rf_tuning/rf_minCR', minCR, '_Allsites.grid.search.rda'))

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)

rf.params <- list()
rf.params$Allsites <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
#rf.params$Allsites$mtry <- 63
#rf.params$Allsites$sampsize <- 35
save(rf.params, file = paste0('R/rf_tuning/rf_optim_params_minCR', minCR, '.rda'))

# RFsites ------------------------------------------------------------------
sites <- readRDS('R/rf_tuning/rf_site_importance_Allsamps.rds') |>
  filter(pval <= 0.05) |>
  pull('loc.site')

param.df  <- rf.param.grid.search(model.df, sites)
save(param.df, file = paste0('R/rf_tuning/rf_minCR', minCR, '_RFsites.grid.search.rda'))

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)

rf.params$RFsites <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
#rf.params$RFsites$mtry <- 12
rf.params$RFsites$sampsize <- 49
save(rf.params, file = paste0('R/rf_tuning/rf_optim_params_minCR', minCR, '.rda'))

# glmnet.5 sites ------------------------------------------------------------------
sites <- readRDS('R/glmnet/glmnet.chosen.sites.rds')$minCR2$alpha.5

param.df  <- rf.param.grid.search(model.df, sites)
save(param.df, file = paste0('R/rf_tuning/rf_minCR', minCR, '_glmnnet.5.grid.search.rda'))

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)

rf.params$glmnet.5 <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
#rf.params$RFsites$mtry <- 6
#rf.params$RFsites$sampsize <- 35
save(rf.params, file = paste0('R/rf_tuning/rf_optim_params_minCR', minCR, '.rda'))

# gamsites ------------------------------------------------------------------
sites <- readRDS('R/gam/gam_significant_sites.rds') |>
  filter(r.sq >= 0.35) |>
  pull('loc.site')

param.df  <- rf.param.grid.search(model.df, sites)
save(param.df, file = paste0('R/rf_tuning/rf_minCR', minCR, '_gamsites.grid.search.rda'))

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)

rf.params$gamsites <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
#rf.params$RFsites$mtry <- 4
#rf.params$RFsites$sampsize <- 35
save(rf.params, file = paste0('R/rf_tuning/rf_optim_params_minCR', minCR, '.rda'))
