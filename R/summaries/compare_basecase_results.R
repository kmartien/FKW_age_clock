library(tidyverse)
library(ggplot2)
library(gridExtra)
source("R/plotting.funcs.R")
load("data/color.palettes.rda")
load('data/age_and_methylation_data.rdata')

minCRs <- c(2, 4)
sites.2.use <- c('RFsites', 'glmnet.5')
weight <- c('CR', 'inv.var', 'sn.wt', 'none')

pred <-
  do.call(rbind, lapply(c('glmnet', 'rf', 'svm', 'gam'), function(m){
    do.call(rbind, lapply(minCRs, function(cr){
      do.call(rbind, lapply(sites.2.use, function(s){
        if (m == 'svm'){
          readRDS(paste0('R/', m, '/', m, '_best_minCR', cr, '_', s, '_ln.rds')) |> 
            mutate(weight = 'none')|> 
            mutate(sites = s)
        } else do.call(rbind, lapply(weight, function(wt){
          readRDS(paste0('R/', m, '/', m, '_best_minCR', cr, '_', s, '_ln_', wt, '.rds')) |> 
            mutate(weight = wt)|> 
            mutate(sites = s)
        })) 
      }))|> 
        select(c('swfsc.id', 'age.resp', 'age.pred', 'resid', 'dev', 'sq.err','weight', 'sites')) |> 
        rename(age.best = age.resp) |>
        mutate(method = m) |> 
        mutate(minCR = cr) 
    }))
  }))

sites.2.use <- c('Allsites')
pred <- bind_rows(pred, 
  do.call(rbind, lapply(c('glmnet', 'rf', 'svm'), function(m){
    do.call(rbind, lapply(minCRs, function(cr){
      do.call(rbind, lapply(sites.2.use, function(s){
        if (m == 'svm'){
          readRDS(paste0('R/', m, '/', m, '_best_minCR', cr, '_', s, '_ln.rds')) |> 
            mutate(weight = 'none')|> 
            mutate(sites = s)
        } else do.call(rbind, lapply(weight, function(wt){
          readRDS(paste0('R/', m, '/', m, '_best_minCR', cr, '_', s, '_ln_', wt, '.rds')) |> 
            mutate(weight = wt)|> 
            mutate(sites = s)
        })) 
      }))|> 
        select(c('swfsc.id', 'age.resp', 'age.pred', 'resid', 'dev', 'sq.err','weight', 'sites')) |> 
        rename(age.best = age.resp) |>
        mutate(method = m) |> 
        mutate(minCR = cr) 
    }))
  }))
) |> 
mutate(model = paste('minCR', minCR, sites, weight, 'best', sep = '_'))

MAE.all <- left_join(
  left_join(
    pred |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(method, sites, weight, minCR) |> 
      summarise(MAE = median(dev)),
    pred |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(method, sites, weight, minCR) |> 
      summarise(lci = quantile(dev, probs = c(.25)))
  ),
  pred |> 
    left_join(age.df) |> 
    filter(age.confidence %in% c(4,5)) |> 
    group_by(method, sites, weight, minCR) |> 
    summarise(uci = quantile(dev, probs = c(.75)))
)
write.csv(MAE.all, file = 'R/summaries/MAE.Allsites.csv')

models.2.plot <- c(
  'minCR_4_RFsites_none_best',
  'minCR_4_Allsites_none_best',
  'minCR_4_glmnet.5_none_best',
  'minCR_2_RFsites_none_best', 
  'minCR_2_RFsites_CR_best',
  'minCR_2_RFsites_sn.wt_best',
  'minCR_2_RFsites_inv.var_best',
  'minCR_4_RFsites_none_RanAge',
  'minCR_4_RFsites_none_RanAgeMeth'
)

pred.2.plot <- filter(pred, model %in% models.2.plot) |> 
  mutate(
    model.name = case_when(
      model == 'minCR_4_RFsites_none_best' ~ 'Base',
      model == 'minCR_4_Allsites_none_best' ~ 'Allsites',
      model == 'minCR_4_glmnet.5_none_best' ~ 'glmnetsites',
      model == 'minCR_2_RFsites_none_best' ~ 'unweighted',
      model == 'minCR_2_RFsites_CR_best' ~ 'CR',
      model == 'minCR_2_RFsites_sn.wt_best' ~ 'sn.wt',
      model == 'minCR_2_RFsites_inv.var_best' ~ 'inv.var',
      model == 'minCR_4_RFsites_none_RanAge' ~ 'RanAge',
      model == 'minCR_4_RFsites_none_RanAgeMeth' ~ 'RanAgeMeth',
      .default = 'other'
    )
  )
box.plot <- 
  pred.2.plot |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = model.name, y = dev)) +
  labs(x = "Model",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
  #  facet_wrap(~method, nrow = 1)
  facet_wrap(~method, nrow = 2,labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/model.boxplot.jpg', width = 1200, height = 1200)
box.plot
dev.off()
# pred <- do.call(rbind, lapply(c('svm', 'rf', 'glmnet', 'gam'), function(m){
#   do.call(rbind, lapply(minCRs, function(cr){
#     readRDS(paste0('R/', m, '/', m, '_best_minCR', cr, '_', sites.2.use, '_ln.rds')) |>
#       select(c('swfsc.id', 'age.resp', 'age.pred', 'resid', 'dev', 'sq.err')) |> 
#       rename(age.best = age.resp) |>
#       mutate(method = m) |> 
#       mutate(minCR = cr)
#   }))
# }))

# Base-case: minCR = 4, weight = none, compare methods and site selection

base.res <-   pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  filter(weight == 'none' & minCR == 4)
  
MAE.base <- 
  base.res |> 
  group_by(method, sites, weight, minCR) |> 
  summarise(MAE = median(dev))

box.plot <- 
  base.res |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = sites, y = dev)) +
  labs(x = "CpG site selection",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
#  facet_wrap(~method, nrow = 1)
  facet_wrap(~method, nrow = 1,labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/site.selection.jpg', width = 900, height = 1200)
box.plot
dev.off()

# sample selection and weighting (sites = RFsites)

minCR.res <-   pred |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  filter(sites == 'RFsites') |> 
  filter(minCR == 2 | weight == 'none' & minCR == 4)

MAE.minCR <- 
  minCR.res |> 
  group_by(method, sites, weight, minCR) |> 
  summarise(MAE = median(dev))

box.plot <- 
  minCR.res |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  mutate(minCR.wt = paste0('minCR', minCR, '_', weight)) |> 
  ggplot() +
  geom_boxplot(aes(x = minCR.wt, y = dev)) +
  labs(x = "minCR",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
#    facet_wrap(~method, nrow = 1)
  facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/minCR.jpg', width = 1800, height = 900)
box.plot
dev.off()

# histogram of MAE by method ------------------------------------------------
MAE.hist <- MAE.all |> 
  ggplot() +
  geom_histogram(aes(x = MAE, fill = method)) 
MAE.hist
-----------------------------------------------------------------------------

plots <- lapply(1:length(dat), function(i){
  p <- plot.loov.res(dat[[i]], min.CR = 4)
  p$p.loov$labels$title <- paste0(names(dat)[i])
  return(p$p.loov)
})

#plots[[1]]$labels$title <- "ENR optimized alpha = 0.1"
#plots[[2]]$labels$title <- "ENR alpha = 0.5"
plots$nrow <- 2
jpeg(file = paste0("R/summaries/plainjane.regression.plots--minCR", minCR, "_", sites.2.use, ".jpg"), width = 1000, height = 800)
do.call(grid.arrange, plots)
dev.off()

# box plots
age.errors.long <- do.call(bind_rows, lapply(1:length(dat), function(i){
  bind_cols(method = names(dat)[i], dat[[i]])
}))
age.errors.long$CR <- as.factor(age.errors.long$age.confidence)
box.plot <- 
  ggplot(age.errors.long) +
  geom_boxplot(aes(x = CR, y = dev)) +
  labs(x = "Confidence Rating",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
  #  facet_wrap(~training.cr, nrow = 3, labeller = labeller(training.cr = training.cr.labs))  
  facet_wrap(~method, nrow = 2)  

# MAE for CR=5 males vs. females
MAE.by.sex <- do.call(rbind, lapply(1:length(dat), function(i){
  filter(dat[[i]], age.confidence == 5) %>% group_by(sex) %>% 
    summarise(MAE = median(dev)) %>% bind_cols(method = names(dat)[i])
}))

# MAE by age and CR
breaks <- c(0,10,25,40)
MAE.by.age <- do.call(cbind, lapply(1:length(dat), function(i){
  do.call(rbind, lapply(1:(length(breaks)-1), function(a){
    filter(dat[[i]], age.best >= breaks[a]) %>% filter(age.best < breaks[a+1]) %>% 
      # group_by(age.confidence) %>% #use this line to compare across CRs
      filter(age.confidence >= 4) %>% #use this line to combined CR4and5
      summarise(MAE = median(dev)) %>% bind_cols(age_bins = breaks[a])
  }))
}))
write.csv(MAE.by.age, file = paste0("R/summaries/MAE.by.age.across.methods-minCR", minCR, "_", sites.2.use, ".csv"))

# compare duplicates
age.df$date.biopsy <- as.POSIXct(age.df$date.biopsy, format = '%Y-%m-%d %H:%M:%S') %>% 
  as.Date()
dupe.sum <- do.call(bind_rows, lapply(1:length(dat), function(m){
  dupes <- left_join(dat[[m]], select(age.df, c(crc.id, swfsc.id, date.biopsy))) %>% 
    filter(n() > 1, .by = crc.id)
  do.call(rbind, lapply(unique(dupes$crc.id), function(i){
    inds <- filter(dupes, crc.id == i) %>% arrange(age.best)
    actual.diff <- difftime(inds$date.biopsy[2], inds$date.biopsy[1], units = "days")/365 %>% 
      as.numeric()
    predicted.diff <- inds$age.pred[2] - inds$age.pred[1]
    return(data.frame(method = names(dat)[m], crc.id = i, actual.diff = actual.diff, predicted.diff = predicted.diff, age.confidence = as.character(inds$age.confidence[1])))
  }))
}))
pair.plot <- ggplot(dupe.sum) +
  geom_point(aes(x = actual.diff, y = predicted.diff, colour = age.confidence), size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = conf.colors) +
  labs(x = "Actual age difference", y = "Predicted age difference") +
  theme(text = element_text(size = 24)) +
  facet_wrap(~method, nrow = 2)
jpeg(filename = paste0("R/summaries/pair.plot_minCR", minCR, "_", sites.2.use, ".jpg"), width = 960, height = 960)
pair.plot
dev.off()

