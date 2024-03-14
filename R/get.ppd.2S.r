#' get.ppd.2S
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
get.ppd.2S <- function(mod, pst.samples=10^3, perc = seq(0, 1, 0.01), type = 'x', ppd.type = 'percentiles', quant = NULL) {
  Vobs = mod$Vobs
  vanilla = mod$vanilla
  max.fu = max( sapply(Vobs, function(x) if(is.infinite(x[length(x)])) return(x[length(x)-1]) else return(x[length(x)]) ))
  if(is.null(quant)) quant = seq(0, max.fu, max.fu/1000 )
  par.list.X = (mod$par.X.bi)
  Z.X = as.matrix(mod$Z.X)
  dist.X = mod$dist.X
  g.fixed = numeric(length = length(Vobs))
  for(i in 1:length(mod$Vobs)) if( length(mod$Vobs[[i]]) == 1 ) g.fixed[i] = 1
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  if(vanilla)  ppd = sample.ppd.vanilla(par.list = par.list.X, Z.X = Z.X, dist.X = dist.X, s = s)
  if(!vanilla) ppd = sample.ppd.xstar(par.list = par.list.X, Z.X = Z.X, Z.W = as.matrix(mod$Z.W),
                                             dist.X = dist.X, s = s, g.fixed = g.fixed, type = type)
  
  degenerate <- apply(ppd, 2, function(x) sum(is.na(x)) == length(x) )
  if(sum(degenerate)>0){
    ppd <- ppd[,!degenerate]
    warning('The prevalence model produced posterior predicted values equal to 1 for all n. This probably signifies a convergence issue.')
  }
  
  ret = list()
  if(ppd.type == 'percentiles'){
    perc = apply( ppd, 2, function(x) {e = ecdf(x); e(quant)} )
    ret$med.cdf = apply(perc, 1, median)
    ret$med.cdf.ci = apply(perc, 1, quantile, c(0.025, 0.975))
    ret$quant = quant
  }
  if(ppd.type == 'quantiles'){
    q = apply( ppd, 2, quantile, perc )
    ret$med.cdf = apply(q, 1, median)
    ret$med.cdf.ci = apply(q, 1, quantile, c(0.025, 0.975))
    ret$perc = perc
  }
  ret
}


