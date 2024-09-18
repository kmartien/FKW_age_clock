rm(list = ls())
library(tidyverse)
library(glmnet)
load("../age_and_methylation_data.rdata")

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
ageCVglmnetLOO <- function(df, sites, alpha) {
  cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df$age.best,
    alpha = alpha,
    foldid = 1:nrow(df)
  ) 
}

# return median cvm at minimum lambda for alpha and nrep replicates
median.cvm.min <- function(alpha, df, nrep) {
  parallel::mclapply(1:nrep, function(i) {
    cv.fit <- tryCatch({
      ageCVglmnetLOO(df, sites.to.keep, alpha)
    }, error = function(e) NULL)
    if(is.null(cv.fit)) NA else {
      cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.min]
    } 
  }, mc.cores = 6) |> 
    unlist() |> 
    median(na.rm = TRUE)
}

# find optimal alpha (lowest median cvm at minimum lambda)
optimum.alpha <- optim(
  par = c(alpha = 0.3),
  fn = median.cvm.min,
  method = "Brent",
  lower = 0.0001,
  upper = 0.5,
  df = filter(model.df, age.confidence %in% 4:5),
  nrep = 1000
)
optimum.alpha

# multiple LOO runs of the glmnet model with alpha at lowest cvm
glmnet.fit <- parallel::mclapply(1:1000, function(i) {
  tryCatch({
    ageCVglmnetLOO(
      filter(model.df, age.confidence %in% 4:5), 
      sites.to.keep, 
      optimum.alpha$par
    )
  }, error = function(e) NULL)
}, mc.cores = 10)

save.image('glmnet.rdata')



# distribution of minimum lambda across model replicates
lambda.min <- sapply(glmnet.fit, function(x) {
  if(is.list(x)) x$lambda.min else NA
})
hist(lambda.min)
swfscMisc::distSmry(lambda.min, method = 'venter')


# distribution of predicted age across model replicates
pred.age <- lapply(glmnet.fit, function(x) {
  predict(
    x, 
    as.matrix(model.df[, sites.to.keep]),
    s = 'lambda.min'
  ) |> 
    as.data.frame() |> 
    rownames_to_column('swfsc.id') |> 
    rename(pred.age = 'lambda.min')
}) |> 
  bind_rows() 

pred.age |> 
  left_join(
    age.df |> 
      select(swfsc.id, age.best, age.confidence),
    by = 'swfsc.id'
  ) |> 
  mutate(age.confidence = as.character(age.confidence)) |> 
  group_by(age.confidence) |> 
  summarize(mse = mean((age.best - pred.age) ^ 2), .groups = 'drop') |> 
  bind_rows(
    pred.age |> 
      left_join(age.df, by = 'swfsc.id') |> 
      summarize(mse = mean((age.best - pred.age) ^ 2)) |> 
      mutate(age.confidence = 'All')
  )


pred.age |>   
  group_by(swfsc.id) |> 
  summarize(
    median.pred.age = median(pred.age),
    lci = quantile(pred.age, 0.025),
    uci = quantile(pred.age, 0.975),
    .groups = 'drop'
  ) |> 
  left_join(
    age.df,
    by = 'swfsc.id'
  ) |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot(aes(x = age.best)) + 
  geom_abline(intercept = 0, slope = 1) +  
  geom_segment(aes(x = age.best, xend = age.best, y = lci, yend = uci, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = median.pred.age, yend = median.pred.age, color = age.confidence), alpha = 0.5) +
  geom_point(aes(y = median.pred.age, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC best age', y = 'GLMNET predicted age') +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(1, 1),
    legend.justification = c(1, 1)
  )
