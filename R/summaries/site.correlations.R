library(tidyverse)
source('R/misc_funcs.R')
load('data/age_and_methylation_data.rdata')

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  )

rf.sites <- selectCpGsites('RFsites')
glmnet.sites <- selectCpGsites('glmnet.5')

age.corr.coeff <- 
  locus.map |> 
  right_join(
  lapply(sites.to.keep, function(s){
    temp <- cor.test(model.df[[s]], model.df$age.best, method = 'pearson', use = 'na.or.complete')
    return(data.frame(loc.site = s, corr.coeff = temp$estimate, p.val = temp$p.value))
  }) |> bind_rows()
) |> 
  mutate(sig = ifelse(p.val <= 0.05, 1, 0))# |> 
  #mutate(sig = as.factor(sig))

jpeg(file = 'R/summaries/site.correlation.with.age.jpg', width = 900, height = 600)
ggplot(age.corr.coeff) +
  geom_point(aes(x = site, y = corr.coeff), size = 1, colour = 'gray50') +
  geom_point(data = filter(age.corr.coeff, p.val <= 0.05), aes(x = site, y = corr.coeff), size = 2, colour = 'black') +
  geom_point(data = filter(age.corr.coeff, loc.site %in% rf.sites), aes(x = site, y = corr.coeff), size = 3, colour = 'red') +
  labs(x = "Site position", y = "Correlation coefficient") +
  theme(text = element_text(size = 15)) +
  facet_wrap(~locus, scales = 'free_x')
dev.off()
