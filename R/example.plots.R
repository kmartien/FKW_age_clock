rm(list=ls())
library(ggplot2)
library(tidyverse)
library(grid)
library(gridExtra)
library(sn)

load("data/age_and_methylation_data.rdata")
age.df <- age.df %>% 
  mutate(label = paste0(crc.id, ", ", swfsc.id)) %>% 
  arrange(crc.id, swfsc.id)

inputIDs <- c("HIPc114","HIPc102","HIPc338")
inds <-  which(age.df$crc.id %in% inputIDs)[-c(1,4)]
#inputIDs <- c("HIPc313","HIPc155")
#inds <-  which(age.df$crc.id %in% inputIDs)
#don't use HIPc173, HIPc190, or HIPc382; data in Eric's dataframe is out of date
#I originally used HIPc313 as the CR=3 example, but then its CR was revised
# to CR=4, so I replaced it with HIPc101, but then it got out of date, so I switched to HIPc102.
p.vec <- c(0.05,0.2,0.55,0.75,1)

plots <- lapply(inds, function(i){
  p <- p.vec[age.df$age.confidence[i]]

  min.age <- age.df$age.min[i]
  max.age <- age.df$age.max[i]
  dp <- unlist(age.df[i, c("sn.location", "sn.scale", "sn.shape")])
  age.vec <- seq(min.age, max.age, length.out = 1000)
  dens.unif <- 1 / (max.age - min.age)
  dens.age <- dsn(age.vec, dp = dp)
  dens.diff <- dens.unif - dens.age

  dens.df <- sapply(c(0,p,1), function(x) dens.unif - (dens.diff * x)) %>% 
    cbind(age.vec) %>% 
    as.data.frame() %>%
    setNames(c("0",age.df$age.confidence[i],"5", "age")) %>% 
    pivot_longer(-age, names_to = "Confidence", values_to = "density") %>% 
    mutate(this.cr = Confidence == as.character(age.df$age.confidence[i])) %>% 
    arrange(Confidence, age) |> 
    #eliminate the next line if you want to plot SN, Uniform, and composite distributions
    filter(this.cr == TRUE)

  y.lims <- c(0,max(dens.age))
  
  line_palette <- c("gray","black","black")
  g <- ggplot(dens.df, aes(age, density)) +
    geom_vline(xintercept = age.df$age.best[i], size = 0.5) +
    geom_vline(xintercept = unlist(age.df[i, c("age.min", "age.max")]), linetype = "solid", size = 0.5) +
    geom_line(aes(
      color = Confidence, 
      linetype = ifelse(Confidence == 5,"solid","dashed") #this seems backward, but it works...
    ), size = 1, show.legend = FALSE) +
    scale_color_manual(values=conf.colors) +
    guides(linetype = "none") +
    labs(x = "Age", y = "Density") + #, title = paste0(age.df$crc.id[i],"\nCR = ",age.df$age.confidence[i])) +
    scale_x_continuous(breaks = unlist(age.df[i,c("age.min", "age.best","age.max")]),
                       labels = as.character(unlist(age.df[i,c("age.min", "age.best","age.max")]))) +
    ylim(y.lims) +
#    theme_minimal() + 
    theme(plot.title =  element_text(size=30), axis.title.x =  element_text(size=25), 
          axis.title.y = element_blank(), axis.text = element_text(size=25))

})
plots <- plots[c(3,1,2)]
plots$nrow <- 3
plots$left <- textGrob("Density", rot = 90, gp = gpar(fontsize = 30))
jpeg(file="example.age.prob.dist.plots.2024FEB.jpg",width = 800, height = 1200)
do.call(grid.arrange, plots)
dev.off()

