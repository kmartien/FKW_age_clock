library(tidyverse)
library(glmnet)
load("data/age_and_methylation_data.rdata")

site.incl.threshold <- 0.8
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
median.cvm.min <- function(alpha, df, nrep) {
  parallel::mclapply(1:nrep, function(i) {
    cv.fit <- tryCatch({
      CVglmnet10K(df, sites.to.keep, alpha)
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

  # find optimal alpha (lowest median cvm at minimum lambda)
  # opt.alpha <- optim(
  #   par = c(alpha = 0.3),
  #   fn = median.cvm.min,
  #   method = "Brent",
  #   lower = 0.0001,
  #   upper = 0.5,
  #   df = train.df,
  #   nrep = 1000
  # )
  label <- paste0('minCR', minCR)
  opt.alpha <- optimum.alpha[[label]]
  
  # multiple LOO runs of the glmnet model with alpha at lowest cvm
  fit <- parallel::mclapply(1:nrep, function(i) {
    tryCatch({
      CVglmnet10K(
        train.df, 
        sites.to.keep, 
        opt.alpha$par
      )
    }, error = function(e) NULL)
  }, mc.cores = 10)
  
  return(list(optimum.alpha = opt.alpha, fit = fit))
})
names(glmnet.optim) <- c('minCR2', 'minCR4')

lapply(glmnet.optim, function(i){i$optimum.alpha}) |>
  saveRDS(file = 'R/glmnet/optimum.alpha.rds')

lapply(glmnet.optim, function(minCR){
  do.call(bind_rows,
  lapply(minCR$fit, function(i){
    coeffs <- coef(i, s = 'lambda.min') 
    data.frame(sites = coeffs@Dimnames[[1]][coeffs@i+1])
  })) |> 
    group_by(sites) |> 
    summarise(count = n()) |> 
    filter(count >= (site.incl.threshold * nrep) & sites != '(Intercept)') |>
    pull(sites)
}) |>
  saveRDS(file = 'R/glmnet/glmnet.chosen.sites.rds')
