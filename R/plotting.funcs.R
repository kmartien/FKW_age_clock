library(ggplot2)
library(gridExtra)
library(viridis)

# Age distribution of training data set
p.age.distribution <- function(dat, min.CR){
  ggplot(filter(dat, age.confidence >= min.CR), aes(x = age.best, fill = as.character(age.confidence))) +
    geom_histogram(bins = 35, binwidth = 1, color = "black") +
    scale_fill_manual(values = conf.palette, name = "Confidence") +
    labs(x = "Age.best", y = "Count") +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.6, 0.8)
    )
}

# Bar chart of social cluster membership by sex
p.cluster.distribution <- function(dat) {
  dat$cluster.louvain <- as.numeric(dat$cluster.louvain)
  ggplot(dat, aes(x = cluster.louvain, col = I("grey"), fill = sex)) +
    geom_histogram(binwidth = 1) +
    scale_fill_manual(values = sex.palette) +
    labs(x = "Social cluster", y = "Count") +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.9, 0.8)
    )
}

# plot loov results by min.CR
plot.loov.res <- function(loov.res, min.CR) {# expected argument is loov.res
  dat <- filter(loov.res, age.confidence >= min.CR)
  med.dev <- dat |> summarise(MAE = round(median(dev),2))
  cor.coeff <- round(cor.test(dat$age.best, dat$age.pred, method = "pearson")$estimate,2)
  loov.regress <- do.call(rbind,lapply(min.CR:5, function(cr){
    dat.filtered <- filter(dat, age.confidence == cr)
    regr <- lm(age.pred~age.best, data = dat.filtered)$coefficients
    names(regr) <- c("intercept", "slope")
    return(regr)
  }))
  fit.sum <- bind_cols(MAE = med.dev, Corr = cor.coeff, loov.regress)
  
  p.loov <- ggplot(dat, aes(x = age.best, y = age.pred)) +
    geom_point(aes(col = as.character(age.confidence)), size = 3) +
    stat_smooth(aes(color = 'black', fill = NULL), 
                show.legend = FALSE,
                method = "lm", 
                formula = y ~ x, 
                geom = "smooth") +
    geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
    labs(x = "Age.best", y = "Predicted age") +
    annotation_custom(tableGrob(fit.sum[1, c('MAE', 'Corr')], theme = ttheme_minimal(base_size = 16), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
    scale_color_manual(values = conf.palette, name = "Confidence") +
    xlim(0,40) + ylim(0,80) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      legend.position = c(0.2, 0.8),
      plot.background = element_rect(color = "black", linewidth = 1)
    )
  # for(i in 1:nrow(fit.sum)){
  #   cr <- i+(min.CR-1)
  #   p.loov <- p.loov + 
  #     geom_abline(slope = fit.sum$slope[i], intercept = fit.sum$intercept[i], color = conf.palette[cr])
  # }
  return(list(p.loov = p.loov, fit.sum = fit.sum))
}

#scatterplot of residuals
plot.residuals <- function(loov.res, min.CR) {# expected argument is loov.res
  dat <- filter(loov.res, age.confidence >= min.CR)
  med.dev <- dat |> summarise(MAE = round(median(dev),2))
  cor.coeff <- round(cor.test(dat$age.best, dat$age.pred, method = "pearson")$estimate,2)
  # loov.regress <- do.call(rbind,lapply(min.CR:5, function(cr){
  #   dat.filtered <- filter(dat, age.confidence == cr)
  #   regr <- lm(age.pred~age.best, data = dat.filtered)$coefficients
  #   names(regr) <- c("intercept", "slope")
  #   return(regr)
  # }))
  fit.sum <- bind_cols(MAE = med.dev, Corr = cor.coeff)
  
  p.loov <- ggplot(dat, aes(x = age.best, y = resid)) +
    geom_point(aes(col = as.character(age.confidence)), size = 3) +
    stat_smooth(aes(color = 'black', fill = NULL), 
                show.legend = FALSE,
                method = "lm", 
                formula = y ~ x, 
                geom = "smooth") +
    geom_abline(slope = 0, color = "black", linewidth = 0.5, linetype = 2) +
    labs(x = "Age.best", y = "Residual") +
    annotation_custom(tableGrob(fit.sum[1, c('MAE', 'Corr')], theme = ttheme_minimal(base_size = 16), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
    scale_color_manual(values = conf.palette, name = "Confidence") +
 #   xlim(0,40) + ylim(0,80) +
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      legend.position = c(0.2, 0.2),
      plot.background = element_rect(color = "black", linewidth = 1)
    )
  return(p.loov)
}

# Distribution of residuals for glmnet
plot.deviation <- function(dat, min.CR){
  dat$age.range <- dat$age.max - dat$age.min
  dat$age.confidence <- as.character(dat$age.confidence)
  ggplot(filter(dat, age.confidence >= min.CR)) +
    geom_point(aes(x = resid, y = age.range, col = age.confidence)) +
    scale_colour_manual(values = conf.palette) +
    labs(x = "Predicted age - CRC age.best", y = "Age.max - age.min") +
    ggtitle(paste0("glmnet Training (CR >=", min.CR, ")")) +
    theme_minimal() +
    theme(
      text = element_text(size = 15)
    )
}

# Histogram of age residuals from LOOCV results
loov.hist <- function(dat, min.cr){ # expected argument is loov.res
  ggplot(filter(dat, age.confidence >= min.cr), aes(x = resid, fill = as.character(age.confidence))) +
    geom_histogram(bins = 35, binwidth = 1, color = "black") +
    scale_fill_manual(values = conf.palette, name = "Confidence") +
    labs(x = "Predicted age - Age.best", y = "Count") +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.2, 0.8)
    )
}

