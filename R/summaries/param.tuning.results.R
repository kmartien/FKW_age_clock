library(tidyverse)
load('R/svm/svm.tuning.rda')
optim.alpha <- readRDS('R/glmnet/optim.alpha.rds')

opt.params <- bind_cols(
lapply(c('minCR2', 'minCR4'), function(minCR){
  lapply(c('Allsites', 'RFsites', 'glmnet.5'), function (s){
    return(data.frame(minCR = minCR, sites = s, alpha = optim.alpha[[minCR]][[s]]))
  }) |> bind_rows()
}) |> bind_rows(),

lapply(c('Allsamps', 'CR4_5'), function(minCR){
  lapply(c('Allsites', 'RFsites', 'glmnet.5'), function (s){
    return(data.frame(minCR = minCR, sites = s, cost = svm.tuning[[minCR]][[s]]$best.parameters$cost, gamma = svm.tuning[[minCR]][[s]]$best.parameters$gamma))
    }) |> bind_rows()
  }) |> bind_rows(),

lapply(c('minCR2', 'minCR4'), function(minCR){
  load(paste0('R/rf_tuning/rf_optim_params_', minCR, '.rda'))
  lapply(c('Allsites', 'RFsites', 'glmnet.5'), function (s){
    return(data.frame(minCR = minCR, sites = s, mtry = rf.params[[s]]$mtry, sampsize = rf.params[[s]]$sampsize))
  }) |> bind_rows()
}) |> bind_rows()
)

write.csv(opt.params, file = 'R/summaries/opt.params.csv')
