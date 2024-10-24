# Plot the results of the base model for each method

df <- pred |> 
  filter(method == 'gam' & sites == 'RFsites' & minCR == 'minCR4' & resample == 'best') |> 
  mutate(weight = factor(weight, levels = c('none', 'CR', 'ci.wt'), ordered = TRUE)) |> 
  left_join(age.df)
base.plot <- lapply(c('none', 'CR', 'ci.wt'), function(w){
  p <- plot.loov.res(filter(df, weight == w), 4)$p.loov + 
    ggtitle(w) 
  return(p)
})
#names(base.plot) <- c('GAM', 'SVM', 'ENR', 'RF')
base.plot$nrow = 4
jpeg(file = 'R/summaries/best.GAM.scatterplot.jpg', width = 600, height = 1200)
do.call(grid.arrange, base.plot)
dev.off()

residuals.plot <-
  filter(df, age.confidence > 3) |> 
  ggplot(aes(x = age.best, y = resid)) +
  geom_point(aes(color = factor(age.confidence)), size = 3, show.legend = FALSE) +
  stat_smooth(aes(color = 'black', fill = NULL), 
              show.legend = FALSE,
              method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  scale_color_manual(values = conf.colors) +
  geom_abline(slope = 0, color = "black", linewidth = 0.5, linetype = 2) +
  labs(x = 'Agebest',
       y = 'Residual') +  
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  facet_wrap(~weight, nrow = 4)
jpeg(file = 'R/summaries/best.GAM.residuals.plot.jpg', width = 600, height = 1200)
residuals.plot
dev.off()

meth.dat.long <- pivot_longer(meth.dat, cols = all_of(sites), names_to = 'site')
ggplot(data = meth.dat.long) +
  geom_boxplot(aes(x = value)) +
  geom_point(aes(x = value, y = 0.1, colour = swfsc.id), data = filter(meth.dat.long, swfsc.id %in% c('z0049052', 'z0175877')))+
  facet_wrap(~site, scales = 'free')
