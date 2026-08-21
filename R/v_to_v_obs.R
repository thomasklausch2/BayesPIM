#' Reduce a full screening series to the series observed under imperfect sensitivity
#'
#' Requires at least one finite screening time after the baseline 0 in every
#' element of `v`, i.e. `length(v[[i]]) >= 3`. `gen_data()` guarantees this by
#' drawing the censoring time as the first post-baseline screening time plus a
#' positive increment. Without one, the prevalent/untested branch below evaluates
#' `1:0`, which counts down rather than skipping, and returns an `NA`-padded
#' series instead of `c(0, Inf)`.
#'
#' The final entry of each `p` is the probability that every remaining test is
#' missed. It is written in closed form as `(1 - kappa)^length(p)` rather than as
#' `1 - sum(p)`: the latter is catastrophic cancellation, losing all precision in
#' that term for long screening series and turning negative once it falls below
#' the rounding error in `sum(p)`, which `rmultinom()` rejects.
#' @noRd
v_to_v_obs = function(v, times, kappa, g, baseline_test){

  kap = kappa
  l = numeric()
  n = length(v)
  for(i in 1:n) l[i] = max(cumsum(times[i] > v[[i]]))
  m = sapply(v, length)
  cens = (m-1) == l

  k = numeric()
  # Unit not tested at baseline
  for(i in 1:n){
  if(baseline_test[i] == 0){
    if(g[i] == 0){
      if(!cens[i]){
        p = kap
        if((m-l-2)[i]>0){
          for(j in 1:(m-l-2)[i] ){
            p[j+1] = kap * (1-kap)^j
          }
        }
        p[length(p) + 1] = (1-kap)^length(p)
        k[i] =sum(rmultinom(1,1,p) * 0:(length(p)-1))
      }
      if(cens[i]){
        k[i] = 0
      }}
    # if(g[i] == 1){
    #   k[i] = 0
    # }
    if(g[i] == 1){
      p = numeric()
      for(j in 1:(m-2)[i] ){
        p[j] = kap * (1-kap)^(j-1)
      }
      p[length(p) + 1] = (1-kap)^length(p)
      k[i]   = sum(rmultinom(1,1,p) * 0:(length(p)-1))
    }
  }
  # Unit tested at baseline
  if(baseline_test[i] == 1){
      if(g[i] == 0){
        if(!cens[i]){
          p = kap
          if((m-l-2)[i]>0){
            for(j in 1:(m-l-2)[i] ){
              p[j+1] = kap * (1-kap)^j
            }
          }
          p[length(p) + 1] = (1-kap)^length(p)
          k[i] =sum(rmultinom(1,1,p) * 0:(length(p)-1))
        }
        if(cens[i]){
          k[i] = 0
        }}
      if(g[i] == 1){
        p = numeric()
        for(j in 1:(m-1)[i] ){
          p[j] = kap * (1-kap)^(j-1)
        }
        p[length(p) + 1] = (1-kap)^length(p)
        k[i]   = sum(rmultinom(1,1,p) * 0:(length(p)-1))
      }
    }
  }
  v_obs = list()
  for(i in 1:length(v)){
    if(baseline_test[i]==0){
      if(g[i] == 0) v_obs[[i]] = v[[i]][1: (l[i]+k[i]+1)]
      if(g[i] == 1) v_obs[[i]] = v[[i]][1: (k[i]+2)]
    }
    if(baseline_test[i]==1){
      if(g[i] == 0) v_obs[[i]] = v[[i]][1: (l[i]+k[i]+1)]
      if(g[i] == 1) v_obs[[i]] = v[[i]][1: (k[i]+1)]
    }
  }
  v_obs
}

