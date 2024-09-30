library(tidyverse)
library(glmnet)
load("data/age_and_methylation_data.rdata")

nrep <- 1000

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') 


# standard inverse variance weighted LOO cv.glmnet run at given alpha
CVglmnet10K <- function(df, sites, alpha) {
  cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df$age.best,
    alpha = alpha,
  ) 
}

# return median cvm at minimum lambda for alpha and nrep replicates
median.cvm.min <- function(alpha, df, sites, nrep) {
  parallel::mclapply(1:nrep, function(i) {
    cv.fit <- tryCatch({
      CVglmnet10K(df, sites, alpha)
    }, error = function(e) NULL)
    if(is.null(cv.fit)) NA else {
      cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.min]
    } 
  }, mc.cores = 6) |> 
    unlist() |> 
    median(na.rm = TRUE)
}

glmnet.optim <- lapply(c(2,4), function(minCR){
  train.df <- filter(model.df, age.confidence >= minCR)
  alpha <- lapply(c('Allsites', 'RFsites', 'glmnet.5', 'gamsites'), function(sites.2.use){
    sites <- sites.to.keep
    if(sites.2.use == 'RFsites'){
      # select important sites from Random Forest tuned with all samples
      sites <- readRDS('R/rf_tuning/rf_site_importance_Allsamps.rds') |>
        filter(pval <= 0.05) |>
        pull('loc.site')
    }
    if(sites.2.use == 'glmnet.5'){
      # select chosen sites from glmnet tuned with all samples, alpha = 0.5
      sites <- readRDS('R/glmnet/glmnet.chosen.sites.rds')$alpha.5$minCR2
    }
    if(sites.2.use == 'gamsites'){
      # select chosen sites from gam.by.site
      sites <- readRDS('R/gam/gam_significant_sites.rds') |>
        filter(r.sq >= 0.35) |>
        pull('loc.site')
    }
    
    #  find optimal alpha (lowest median cvm at minimum lambda)
    opt.alpha <- optim(
      par = c(alpha = 0.3),
      fn = median.cvm.min,
      method = "Brent",
      lower = 0.0001,
      upper = 0.5,
      df = train.df,
      sites = sites,
      nrep = nrep
    )
    return(opt.alpha$par)
  })
  names(alpha) <- c('Allsites', 'RFsites', 'glmnet.5', 'gamsites')
  return(alpha)
})
names(glmnet.optim) <- c('minCR2', 'minCR4')
saveRDS(glmnet.optim, file = 'R/glmnet/optim.alpha.rds')

