library(tidyverse)
library(ggplot2)
library(gridExtra)
source("R/plotting.funcs.R")
#load("data/model.params.rda")
load("data/color.palettes.rda")
load('data/age_and_methylation_data.rdata')

#age.transform <- "ln"
minCR <- 4
sites.2.use <- 'RFsites'

dat <- lapply(c('svm', 'rf', 'glmnet'), function(method){
  readRDS(paste0('R/', method, '/', method, '_best_minCR', minCR, '_', sites.2.use, '.rds')) |>
    rename(age.best = age.resp) |>
    left_join(select(age.df, c(swfsc.id, age.confidence, sex)))
})
names(dat) <- c('svm', 'rf', 'glmnet')
dat$gam <- readRDS('R/gam/gam_best_CR4and5.rds') |>
  rename(age.best = age.resp) |>
  left_join(select(age.df, c(swfsc.id, age.confidence, sex)))

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

