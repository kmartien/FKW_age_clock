library(tidyverse)
load('data/age_and_methylation_data.rdata')

sites.2.use <- 'RFsites'
minCR <- 4

method <- 'gam'
res.path <- file.path('R', method, method)

pred <- bind_rows(
  readRDS(paste0(res.path, '_best_CR4and5.rds')) |>
    mutate(type = 'best'),
  readRDS(paste0(res.path, '_ran_age_CR4and5.rds')) |>
    mutate(type = 'ran.age'),
  readRDS(paste0(res.path, '_ran_age_meth_CR4and5.rds')) |>
    mutate(type = 'ran.age.meth')
)
# pred <- bind_rows(
#   readRDS(paste0(res.path, '_best_minCR', minCR, '_', sites.2.use, '.rds')) |>
#     mutate(type = 'best'),
#   readRDS(paste0(res.path, '_ran_age_minCR', minCR, '_', sites.2.use, '.rds')) |>
#     mutate(type = 'ran.age'),
#   readRDS(paste0(res.path, '_ran_age_meth_minCR', minCR, '_', sites.2.use, '.rds')) |>
#     mutate(type = 'ran.age.meth')
# )

errSmry <- function(df) {
  df |> 
    summarize(
      mean.dev = mean(dev),
      median.dev = median(dev),
      mode.dev = modeest::venter(dev),
      mean.se = mean(sq.err),
      median.se = median(sq.err),
      mode.se = modeest::venter(sq.err), 
      cor.age = cor(age.pred, age.resp),
      pct.gt.best = mean(age.pred > age.best),
      .groups = 'drop'
    )
}

smry <- pred |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = as.character(age.confidence)) |> 
  group_by(type, age.confidence) |> 
  errSmry() |> 
  bind_rows(
    pred |> 
      left_join(age.df, by = 'swfsc.id') |> 
      group_by(type) |> 
      errSmry() |> 
      mutate(age.confidence = 'All')
  ) |> 
  arrange(type, age.confidence)
print(smry)
write.csv(smry, file = paste0('R/summaries/', method, '_minCR', minCR, '_', sites.2.use, '.csv'))

# Squared error distribution ----------------------------------------------

pdf(file = paste0('R/summaries/', method, '_minCR', minCR, '_', sites.2.use, '.pdf'))
pred |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() +
  geom_histogram(aes(sq.err, fill = age.confidence)) +
  scale_fill_manual(values = conf.colors) +
  scale_x_log10() +
  facet_grid(type ~ age.confidence, scales = 'free') + 
  theme(legend.position = 'top')


# Distribution of predictions by ID ---------------------------------------

pred |> 
  filter(type == 'best') |> 
  ggplot() + 
  geom_histogram(aes(x = age.pred)) +
  labs(x = paste0(method, ' predicted age'), title = 'best')

pred |> 
  filter(type == 'ran.age') |>
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_histogram(aes(x = age.pred, fill = age.confidence)) +
  geom_vline(
    aes(xintercept = median), 
    data = pred |> 
      filter(type == 'ran.age') |> 
      group_by(swfsc.id) |> 
      summarize(median = median(age.pred), .groups = 'drop') 
  ) +
  facet_wrap(~ swfsc.id, scales = 'free_y') +
  scale_fill_manual(values = conf.colors) +
  labs(x = paste0(method, ' predicted age'), title = 'ran.age') +
  theme(legend.position = 'none')

pred |> 
  filter(type == 'ran.age.meth') |> 
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_histogram(aes(x = age.pred, fill = age.confidence)) + 
  geom_vline(
    aes(xintercept = median), 
    data = pred |> 
      filter(type == 'ran.age.meth') |> 
      group_by(swfsc.id) |> 
      summarize(median = median(age.pred), .groups = 'drop') 
  ) +
  facet_wrap(~ swfsc.id, scales = 'free_y') +
  scale_fill_manual(values = conf.colors) +
  labs(x = paste0(method, ' predicted age'), title = 'ran.age.meth') +
  theme(legend.position = 'none')


# Distribution of residuals by ID -----------------------------------------

pred |> 
  filter(type == 'best') |> 
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_histogram(aes(x = resid, fill = age.confidence)) +
  scale_fill_manual(values = conf.colors) +
  labs(x = 'Residuals', title = 'best')

pred |> 
  filter(type == 'ran.age') |>
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_histogram(aes(x = resid, fill = age.confidence)) +
  geom_vline(
    aes(xintercept = median), 
    data = pred |> 
      filter(type == 'ran.age') |> 
      group_by(swfsc.id) |> 
      summarize(median = median(resid), .groups = 'drop') 
  ) +
  facet_wrap(~ swfsc.id, scales = 'free_y') +
  scale_fill_manual(values = conf.colors) +
  labs(x = 'Residuals', title = 'ran.age') +
  theme(legend.position = 'none')

pred |> 
  filter(type == 'ran.age.meth') |> 
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_histogram(aes(x = resid, fill = age.confidence)) + 
  geom_vline(
    aes(xintercept = median), 
    data = pred |> 
      filter(type == 'ran.age.meth') |> 
      group_by(swfsc.id) |> 
      summarize(median = median(resid), .groups = 'drop') 
  ) +
  facet_wrap(~ swfsc.id, scales = 'free_y') +
  scale_fill_manual(values = conf.colors) +
  labs(x = 'Residuals', title = 'ran.age.meth') +
  theme(legend.position = 'none')


# Predicted vs best age ---------------------------------------------------

pred |> 
  group_by(type, swfsc.id) |> 
  summarize(
    as.data.frame(rbind(HDInterval::hdi(age.pred))), 
    age.pred = median(age.pred),
    .groups = 'drop'
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_abline(intercept = 0, slope = 1) +  
  geom_segment(aes(x = age.best, xend = age.best, y = lower, yend = upper, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = age.pred, yend = age.pred, color = age.confidence), alpha = 0.5) +
  geom_point(aes(x = age.best, y = age.pred, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC age', y = paste0(method, ' predicted age')) +
  facet_grid(age.confidence ~ type) +
  theme(legend.position = 'none')


# Residual vs best age ----------------------------------------------------

pred |> 
  group_by(type, swfsc.id) |> 
  summarize( 
    age.pred = median(age.pred),
    resid.lci = unname(quantile(resid, 0.025)),
    resid.uci = unname(quantile(resid, 0.975)),
    resid = median(resid)
  ) |> 
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_hline(yintercept = 0) + 
  geom_segment(aes(x = age.best, xend = age.best, y = resid.lci, yend = resid.uci, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = resid, yend = resid, color = age.confidence), alpha = 0.5) +
  geom_point(aes(x = age.best, y = resid, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC age', y = 'Residual') +
  facet_grid(age.confidence ~ type, scales = 'free_y') +
  theme(legend.position = 'none')

dev.off()
