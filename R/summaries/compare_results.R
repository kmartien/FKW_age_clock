library(tidyverse)
library(ggplot2)
library(gridExtra)
source("R/plotting.funcs.R")
load("data/color.palettes.rda")
load('data/age_and_methylation_data.rdata')

# subset all data
age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

minCRs <- c(2, 4)
sites.2.use <- c('RFsites', 'glmnet.5', 'Allsites')
wts.2.use <- c('CR', 'ci.wt', 'none')

res.files <- bind_rows(
  lapply(c('svm', 'gam', 'glmnet', 'rf'), function(m){
  fnames <- list.files(paste0('R/', m), pattern = paste0(m, '_best'))
      bind_cols(fnames, do.call(rbind, strsplit(fnames, split = '_')))
}) |> bind_rows(),
lapply(c('svm', 'gam', 'glmnet', 'rf'), function(m){
  fnames <- list.files(paste0('R/', m), pattern = paste0(m, '_ranAge'))
  bind_cols(fnames, do.call(rbind, strsplit(fnames, split = '_')))
}) |> bind_rows()
)
names(res.files) <- c('fname', 'method', 'resample', 'minCR', 'sites', 'age.transform', 'weight')
res.files$weight <- substr(res.files$weight, 1, nchar(res.files$weight) - 4)

pred <- bind_rows(
  left_join(
  filter(res.files, resample == 'best'),
  lapply(1:nrow(res.files), function(f){
    data.frame(fname = res.files$fname[f], readRDS(paste0('R/', res.files$method[f], '/', res.files$fname[f])))
  }) |> bind_rows(),
  by = 'fname'
) |> select(-c(lci, uci, k)) |> 
  filter(weight %in% wts.2.use, sites %in% sites.2.use) |> 
  left_join(select(age.df, c(swfsc.id, age.best))) |> 
  # correcting resid, dev, and sq.err for models with resampling
  mutate(resid = age.pred - age.best,
         dev = abs(resid),
         sq.err = resid^2,
         model = paste(resample, minCR, sites, age.transform, weight, sep = '_')), 

left_join(
  filter(res.files, resample != 'best'),
  lapply(1:nrow(res.files), function(f){
    data.frame(fname = res.files$fname[f], readRDS(paste0('R/', res.files$method[f], '/', res.files$fname[f])))
  }) |> bind_rows(),
  by = 'fname'
) |> select(-c(lci, uci, k)) |> 
  filter(weight %in% wts.2.use, sites %in% sites.2.use) |> 
  group_by(method, age.transform, sites, weight, minCR, resample, swfsc.id) |> 
  summarise(age.pred = quantile(age.pred, probs = c(0.5))) |> 
  left_join(select(age.df, c(swfsc.id, age.best))) |> 
  mutate(resid = age.pred - age.best,
         dev = abs(resid),
         sq.err = resid ^ 2,
         model = paste(resample, minCR, sites, age.transform, weight, sep = '_'))
)

MAE.all <- 
    pred |> 
      left_join(age.df) |> 
      filter(age.confidence %in% c(4,5)) |> 
      group_by(method, sites, weight, minCR, resample) |> 
      summarise(MAE = round(median(dev), 2),
                lci = round(quantile(dev, probs = c(.25)),2),
                uci = round(quantile(dev, probs = c(.75)),2),
                Corr = round(cor.test(age.best, age.pred, method = "pearson")$estimate,2))
write.csv(MAE.all, file = 'R/summaries/MAE.all.csv')

models.2.plot <- c(
  'best_minCR4_RFsites_ln_none',
  'best_minCR4_Allsites_ln_none',
  'best_minCR4_glmnet.5_ln_none',
  'best_minCR2_RFsites_ln_none', 
  'best_minCR2_RFsites_ln_CR',
  'best_minCR2_RFsites_ln_sn.wt',
  'best_minCR2_RFsites_ln_inv.var',
  'ranAge_minCR4_RFsites_ln_none',
  'ranAgeMeth_minCR4_RFsites_ln_none'
)

pred.2.plot <- filter(pred, model %in% models.2.plot) 
pred.2.plot$model <- factor(pred.2.plot$model, levels = models.2.plot, ordered = TRUE) 

# Base model box plot -----------------------------------------------

df <- pred |> 
  filter(model == 'best_minCR4_RFsites_ln_none') |> 
  left_join(age.df)
base.plot <- lapply(c('gam', 'glmnet', 'rf', 'svm'), function(m){
  p <- plot.loov.res(filter(df, method == m), 4)$p.loov + 
    ggtitle(m) 
  #    p$theme$title <- element_text(m)
  return(p)
})
names(base.plot) <- c('svm', 'gam', 'glmnet', 'rf')
base.plot$nrow = 2
jpeg(file = 'R/summaries/base.model.plot.jpg', width = 1200, height = 1200)
do.call(grid.arrange, base.plot)
dev.off()

residuals.plot <- 
  ggplot(filter(df, age.confidence > 3)) +
  # geom_histogram(aes(x = resid, colour = 'gray', fill = factor(age.confidence))) +
  # scale_fill_manual(values = conf.colors) +
  geom_point(aes(x = age.best, y = resid, color = factor(age.confidence))) +
  scale_color_manual(values = conf.colors) +
  geom_abline(slope = 0, intercept = 0) +
  facet_wrap(~method)
residuals.plot

box.plot <- 
  pred.2.plot |> 
  filter(model == 'best_minCR4_RFsites_ln_none') |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  ggplot() +
  geom_boxplot(aes(x = method, y = dev)) +
  scale_x_discrete(labels = c('GAM', 'ENR', 'RF', 'SVM')) +
  labs(x = "Method",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24))
jpeg(file = 'R/summaries/base.model.boxplot.jpg', width = 700, height = 1200)
box.plot
dev.off()

# Site selection box plot -----------------------------------------------

box.plot <- 
  pred.2.plot |> 
  filter(model %in% c('best_minCR4_RFsites_ln_none', 
                      'best_minCR4_Allsites_ln_none',
                      'best_minCR4_glmnet.5_ln_none')) |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  mutate(sites = factor(sites, levels = c('RFsites', 'glmnet.5', 'Allsites'), ordered = TRUE)) |> 
  #  mutate(age.conficence = as.factor(age.confidence)) |> 
  #  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = sites, y = dev)) +
  scale_x_discrete(labels = c('RF sites', 'ENR sites', 'All sites')) +
  #  ylim(c(0,20)) +
  labs(x = "CpG set used",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/site.selection.boxplot.jpg', width = 700, height = 1200)
box.plot
dev.off()


# minCR and weight box plot -----------------------------------------------

box.plot <- 
  pred.2.plot |> 
  filter(model %in% c('best_minCR4_RFsites_ln_none', 
                      'best_minCR2_RFsites_ln_none', 
                      'best_minCR2_RFsites_ln_CR',
                      'best_minCR2_RFsites_ln_sn.wt')) |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  mutate(weight = factor(weight, levels = c('none', 'CR', 'sn.wt'), ordered = TRUE)) |> 
  #  mutate(age.conficence = as.factor(age.confidence)) |> 
  #  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = dev)) +
  scale_x_discrete(labels = c('Base', 'Unweighted', 'CR', 'SN weight')) +
  #  ylim(c(0,20)) +
  labs(x = "Weight scheme",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/minCR.weight.boxplot.jpg', width = 700, height = 900)
box.plot
dev.off()

# summarise RanAge and RanAgeMeth -----------------------------------------------

box.plot <- 
  pred |> 
  filter(model %in% c('best_minCR4_RFsites_ln_none',
                      'ranAge_minCR4_RFsites_ln_none',
                      'ranAgeMeth_minCR4_RFsites_ln_none',
                      'best_minCR2_RFsites_ln_none',
                      'ranAge_minCR2_RFsites_ln_none',
                      'ranAgeMeth_minCR2_RFsites_ln_none')) |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  #  mutate(age.conficence = as.factor(age.confidence)) |> 
  #  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = dev)) +
#  scale_x_discrete(labels = c('HCsamps, best', 'HCsamps, Age', 'Allsamps, Age', 'Allsamps, Age')) +
#  ylim(c(0,25)) +
  labs(x = "Resampling",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~method, nrow = 1, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/resampling.boxplot.jpg', width = 700, height = 900)
box.plot
dev.off()

pred.ran <- 
  left_join(
  filter(res.files, resample != 'best'),
  lapply(1:nrow(res.files), function(f){
    data.frame(fname = res.files$fname[f], readRDS(paste0('R/', res.files$method[f], '/', res.files$fname[f])))
  }) |> bind_rows(),
  by = 'fname'
) |> select(-c(lci, uci, k)) |>
  filter(weight %in% wts.2.use, sites %in% sites.2.use) |>
  group_by(method, age.transform, sites, weight, minCR, resample, swfsc.id) |> 
  summarise(med.pred = quantile(age.pred, probs = c(0.5)),
            lci = quantile(age.pred, probs = c(0.025)),
            uci = quantile(age.pred, probs = c(0.975))) |> 
  left_join(select(age.df, c(swfsc.id, age.best, age.confidence))) |> 
  filter(age.confidence > 3) |> 
  mutate(in.ci = age.best >= lci & age.best <= uci,
         model = paste(resample, minCR, sites, age.transform, weight, sep = '_')) |>
  group_by(method, age.transform, sites, weight, minCR, resample) |>
  summarise(num.in.ci = sum(in.ci))

# compare duplicates -----------------------------------------------

age.df$date.biopsy <- as.POSIXct(age.df$date.biopsy, format = '%Y-%m-%d %H:%M:%S') %>% 
  as.Date()

pair.sum <- lapply(
  select(age.df, c(crc.id, swfsc.id, date.biopsy)) |> 
    filter(n() > 1, .by = crc.id) |> 
    pull(crc.id) |> 
    unique(), 
  function(i){
    inds <- filter(age.df, crc.id == i) %>% arrange(date.biopsy)
    c(
      pair.id = i,
      old.ind = inds$swfsc.id[2], 
      young.ind = inds$swfsc.id[1], 
      age.diff = difftime(inds$date.biopsy[2], inds$date.biopsy[1], units = "days")/365
    )
  }) |> bind_rows() 

age.diff <- lapply(c('gam','glmnet', 'rf', 'svm'), function(mthd){
  mthd.df <- filter(pred, method == mthd & resample == 'best')
  lapply(unique(mthd.df$model), function(mdl){
    df <- filter(mthd.df, model == mdl)
    do.call(rbind, lapply(1:nrow(pair.sum), function(i){
      return(data.frame(
        pair.id = pair.sum$pair.id[i],
        method = mthd,
        model = mdl,
        pred.diff = filter(df, swfsc.id == pair.sum$old.ind[i])$age.pred - filter(df, swfsc.id == pair.sum$young.ind[i])$age.pred
      ))
    }) )
  }) |> bind_rows()
}) |> bind_rows()

pair.plot <- filter(age.diff, model == 'best_minCR4_RFsites_ln_none') |> 
  left_join(pair.sum) |> 
  ggplot() +
  geom_point(aes(x = as.numeric(age.diff), y = pred.diff)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~method, nrow = 2, labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
pair.plot

# MAE by age and CR -----------------------------------------------
breaks <- c(0,10,25,40)
df <- pred |>
  left_join(age.df) |> 
  filter(model == 'best_minCR4_RFsites_ln_none' & age.confidence > 3)

MAE.by.age <-  do.call(rbind, lapply(1:(length(breaks)-1), function(a){
  filter(df, age.best >= breaks[a]) %>% 
    filter(age.best
           < breaks[a+1]) %>% 
    group_by(method) |> 
    summarise(MAE = median(dev)) |> 
    mutate(age_bin = breaks[a]) #%>% 
})) |> 
  pivot_wider(names_from = method, values_from = MAE)
write.csv(MAE.by.age, file = paste0("R/summaries/MAE.by.age.base.models.csv"))

# Check base GAM results ------------------------------------------------

base.gam <- readRDS(file = 'R/gam/gam_best_minCR4_RFsites_ln_none.rds') |> 
  left_join(age.df) |> 
  mutate(in.ci = ifelse(age.best >= lci & age.best <= uci, 1, 0))

filter(base.gam, age.confidence > 3) |> 
  group_by(age.confidence) |> 
  summarise(sum(in.ci))

###########################################################################
# code I'm not currently using 


box.plot <- 
  pred.2.plot |> 
  left_join(age.df) |> 
  filter(age.confidence %in% c(4,5)) |> 
  mutate(age.conficence = as.factor(age.confidence)) |> 
  mutate(minCR = as.factor(minCR)) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = dev)) +
  labs(x = "",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_blank()) +
  facet_wrap(~method, nrow = 2,labeller = labeller(method = c(
    'svm' = 'SVM', 'glmnet' = 'ENR', 'gam' = 'GAM', 'rf' = 'RF'
  )))
jpeg(file = 'R/summaries/model.boxplot.jpg', width = 1200, height = 1800)
box.plot
dev.off()


# plot age.best vs. inv.var
plots <-
  lapply(c('CR', 'sn.wt', 'inv.var'), function(weight){
    age.df |> 
      mutate(wt = if (weight == 'CR') age.confidence else {if(weight == 'inv.var') 1/age.var else confidence.wt}) |> 
#      mutate(age.confidence = factor(age.confidence)) |> 
      ggplot() +
      geom_point(aes(x = age.best, y = wt))#, colour = age.confidence)) +
#      scale_fill_manual(values = conf.colors)
  })



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

#compare duplicates  
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

