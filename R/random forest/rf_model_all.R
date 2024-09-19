rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)

load("data/age_and_methylation_data.rdata")

#' run randomForest over sampsize and mtry grid 
#' and report deviance stats
rf.param.grid.search <- function(model.df, sites){
  # grid of sampsize and mtry to optimize MSE over
  params <- expand.grid(
    sampsize = seq(round(nrow(model.df)*.3), round(nrow(model.df)*.7), by = 1), 
    mtry = seq(round(length(sites)*.1), round(length(sites)*.5), by = 1), 
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
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') |> 
  mutate(inv.var = 1 / age.var)

param.df  <- rf.param.grid.search(model.df, sites.to.keep)
save(param.df, file = "R/random forest/RF.Allsamps.Allsites.grid.search.rda")

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)

rf.params <- list()
rf.params$Allsites <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
#rf.params$Allsites$mtry <- 60
#rf.params$Allsites$sampsize <- 35
save(rf.params, file = 'R/random forest/rf_optim_params_Allsamps.rda')

# run rfPermute at sampsize and mtry with minimum MSE
rp <- rfPermute(
  y = model.df$age.best,
  x = model.df[, sites.to.keep],
  mtry = rf.params$Allsites$mtry,
  ntree = 20000,
  sampsize = rf.params$Allsites$sampsize,
  importance = TRUE,
  replace = FALSE,
  num.rep = 1000,
  num.cores = 10
)

save.image('R/random forest/rf_model_all.rdata')

# save site importance scores and p-values
rp |> 
  importance() |> 
  as.data.frame() |> 
  rownames_to_column('loc.site') |> 
  rename(
    incMSE = '%IncMSE',
    pval = '%IncMSE.pval'
  ) |> 
  select(loc.site, incMSE, pval) |> 
  saveRDS('R/random forest/rf_site_importance_Allsamps.rds')

# select important sites from Random Forest
sites <- readRDS('R/random forest/rf_site_importance_Allsamps.rds') |> 
  filter(pval <= 0.1) |> 
  pull('loc.site')

param.df  <- rf.param.grid.search(model.df, sites)
save(param.df, file = "R/random forest/RF.Allsamps.RFsites.grid.search.rda")

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)

rf.params$RFsites <- filter(param.df, mse == min(param.df$mse)) |>
  select(c(mtry, sampsize))
#rf.params$RFsites$mtry <- 19
#rf.params$RFsites$sampsize <- 62
save(rf.params, file = 'R/random forest/rf_optim_params_Allsamps.rda')

# model diagnostics
print(rp)
plotTrace(rp)
plotImportance(rp, imp = '%IncMSE', alpha = 0.2)
plotImportance(rp, sig = TRUE)

model.df |> 
  mutate(
    age.est = rp$rf$predicted,
    age.confidence = factor(age.confidence)
  ) |> 
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(age.best, age.est, color = age.confidence), size = 3) +
  geom_segment(aes(x = age.min, xend = age.max, y = age.est, yend = age.est, color = age.confidence), alpha = 0.5) +
  scale_color_manual(values = conf.colors) +
  theme(legend.position = 'top')

age.pred <- rp$rf$predicted |> 
  enframe(name = 'swfsc.id', value = 'age.pred') |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(mse = (age.best - age.pred) ^ 2)

age.pred |> 
  mutate(age.confidence = as.character(age.confidence)) |> 
  group_by(age.confidence) |> 
  summarize(mean.mse = mean(mse), .groups = 'drop') |> 
  bind_rows(
    data.frame(age.confidence = 'All', mean.mse = mean(age.pred$mse))
  )