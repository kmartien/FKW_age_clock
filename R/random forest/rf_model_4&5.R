rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)

load("R/data/age_and_methylation_data.rdata")

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep & age.confidence >= 4) |>
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') |> 
  mutate(var.wt = 1 / age.var)

params <- expand.grid(
  sampsize = seq(floor(nrow(model.df)/4), round(nrow(model.df)*0.632), by = 2), 
  mtry = seq(30, 70, by = 2), 
  KEEP.OUT.ATTRS = FALSE
)

#' run randomForest over sampsize and mtry grid 
#' and report deviance stats
param.df <- parallel::mclapply(1:nrow(params), function(i){
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

plot <- param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'YlOrRd')

rf.params <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
save(rf.params, file = 'R/random forest/rf_optim_params_CR4&5_All.rda')

rp <- rfPermute(
  y = model.df$age.best,
  x = model.df[, sites.to.keep],
  mtry = rf.params$mtry,
  ntree = 10000,
  sampsize = rf.params$sampsize,
  importance = TRUE,
  replace = FALSE,
  num.rep = 1000,
  num.cores = 10
)

# save site importance scores and p-values
rf_site_importance_CR4_5 <- 
  rp |> 
  importance() |> 
  as.data.frame() |> 
  rownames_to_column('loc.site') |> 
  rename(
    incMSE = '%IncMSE',
    pval = '%IncMSE.pval'
  ) |> 
  select(loc.site, incMSE, pval)  
  saveRDS(rf_site_importance_CR4_5, 'R/random forest/rf_site_importance_CR4&5.rds')

print(rp)
plotTrace(rp)
plotInbag(rp)
plotImportance(rp)
plotImportance(rp, sig = TRUE)

#' re-run randomForest over sampsize and mtry grid with only important sites
#' and report deviance stats

sites <- rf_site_importance_CR4_5 |> 
  filter(pval <= 0.1) |> 
  pull('loc.site')

params <- expand.grid(
  sampsize = seq(floor(nrow(model.df)/4), round(nrow(model.df)*0.632), by = 2), 
  mtry = seq(floor(length(sites)*0.15), round(length(sites)*.4), by = 1), 
  KEEP.OUT.ATTRS = FALSE
)

#' run randomForest over sampsize and mtry grid 
#' and report deviance stats
param.df <- parallel::mclapply(1:nrow(params), function(i){
  rf <- randomForest(
    y = model.df$age.best,
    x = model.df[, sites],
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

param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'YlOrRd')

rf.params <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
save(rf.params, file = 'R/random forest/rf_optim_params_CR4&5_RFsites.rda')
