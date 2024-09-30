library(tidyverse)
library(randomForest)
library(rfPermute)

load("data/age_and_methylation_data.rdata")

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

#------------------------------------------------------------------------
# All samples
load('R/rf_tuning/rf_optim_params_Allsamps.rda')
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

save.image('R/random forest/rf_model_CR4&5.rdata')

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

#------------------------------------------------------------------------
# CR 4 and 5 samples
load('R/rf_tuning/rf_optim_params_CR4&5.rda')
model.df <- filter(model.df, age.confidence >= 4)
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

save.image('R/random forest/rf_model_CR4&5.rdata')

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
  saveRDS('R/random forest/rf_site_importance_CR4&5.rds')
