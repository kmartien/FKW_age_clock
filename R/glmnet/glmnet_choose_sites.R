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


# standard LOO cv.glmnet run at given alpha
CVglmnet10K <- function(df, sites, alpha) {
  cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df$age.best,
    alpha = alpha,
  ) 
}

glmnet.chosen.sites <- 
  lapply(c(2,4), function(minCR){
    train.df <- filter(model.df, age.confidence >= minCR)
    optimum.alpha <- readRDS('R/glmnet/optimum.alpha.rds')[[paste0('minCR', minCR)]]
    sites <- lapply(c(optimum.alpha$par, 0.5), function(alpha){
      
      # multiple LOO runs of the glmnet model at optimum alpha
      fit <- parallel::mclapply(1:nrep, function(i) {
        tryCatch({
          CVglmnet10K(
            train.df, 
            sites.to.keep, 
            alpha
          )
        }, error = function(e) NULL)
      }, mc.cores = 10)
      
      do.call(bind_rows,
              lapply(fit, function(i){
                coeffs <- coef(i, s = 'lambda.min') 
                data.frame(sites = coeffs@Dimnames[[1]][coeffs@i+1])
              })) |> 
        group_by(sites) |> 
        summarise(count = n()) |> 
        filter(count >= (site.incl.threshold * nrep) & sites != '(Intercept)') |>
        pull(sites)
    })
    names(sites) <- c('alpha.opt', 'alpha.5')
    return(sites)
  })
names(glmnet.chosen.sites) <- c('minCR2', 'minCR4')

saveRDS(glmnet.chosen.sites, file = 'R/glmnet/glmnet.chosen.sites.rds')

