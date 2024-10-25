select.glmnet.chosen.sites <- function(incl.prob.threshold){
  load(file = paste0("data/glmnet_stability_results-sites9_logit_mincov100.rda"))
  filter(site.incl.counts, prob >= incl.prob.threshold)
}