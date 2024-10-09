calc.ci.wt <- function(model.df, logit.meth.normal.params, ndraws = 1000){
  left_join(
    model.df |> 
      mutate(age.confidence = factor(age.confidence)), 
    
    bind_rows(
      lapply(1:ndraws, function(i){
        sampleAgeMeth(model.df, logit.meth.normal.params) |>
          select(swfsc.id, age.ran)
      })) |> 
      group_by(swfsc.id) |> 
      summarise(lci = quantile(age.ran, probs = c(0.25)),
                uci = quantile(age.ran, probs = c(0.75)),
                ci.range = uci - lci) |> 
      select(swfsc.id, ci.range)
  ) |> 
    mutate(ci.wt = 1 - (ci.range / age.range),
           ci.wt = (ci.wt - 0.3) / 0.4) |> 
    pull(ci.wt)
}

