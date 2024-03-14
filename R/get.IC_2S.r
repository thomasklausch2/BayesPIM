#' get.IC_2S
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
get.IC_2S = function(mod, samples = nrow(mod$par.X.bi[[1]]), cores = NULL){
  update.kappa = mod$update.kappa
  vanilla = mod$vanilla
  p1.X = ncol(mod$Z.X)+1
  if(!vanilla) p1.W = ncol(mod$Z.W)
  
  if(is.null(cores)) cores = detectCores()
  m.X = trim.mcmc(mod$par.X.all, burnin = round(nrow(as.matrix(mod$par.X.all[1]))/2)+1)
  m.X = as.matrix(m.X)
  if(!update.kappa) m.X = cbind(m.X, mod$kappa)
  if(samples >  nrow(m.X)) {stop ('more samples than mcmc draws selected')}
  
  ncol.X= ncol(m.X)
  m.X[,(p1.X + 1)] = log(m.X[,(p1.X + 1)])
  m.s = m.X[sample(1:nrow(m.X), samples, replace=F),]
  
  cl    = makePSOCKcluster(cores)
  clusterSetRNGStream(cl)
  registerDoParallel(cl)
  s = round(seq (1, samples, length.out = cores+1))
  pst.mean = apply(m.X, 2, mean) 
  
  run = foreach(j = 1:cores # ndraws=rep(mc, cores)+1:cores
                , .export = c('pdist', 'qdist', 'rdist', 'ddist', 'Lobs_2S',
                              'ploglog', 'rloglog','qloglog', 'dloglog','P_vobs', 'geom', 'geom.inf', 'P_vobs2',
                              'pllogis', 'rllogis', 'dllogis', 'qllogis', 'trans.par', 'trans.par.ind.norm')
                
  ) %dopar% {
    m.s.cores = m.s[s[j]:(s[j+1]-1),]
    out = apply(m.s.cores,1, function(x) Lobs_2S(est=x, mod=mod, log.scale = F, sumup=F) ) 
  }
  stopCluster(cl)
  
  run = do.call(cbind, run)
  
  lppd = sum (log(apply(run,1, mean)))
  q1   = sum (apply( log(run), 1, mean))
  q2   = sum( apply( log(run), 1, var) )
  q3   = mean( apply( log(run),2, sum) )
  q4   = sum(Lobs_2S(pst.mean, mod=mod, log.scale = T, sumup=F))
  DIC  = -2* ( q4 - 2*(q4-q3) )
  WAIC1    = -2*(-lppd + 2*q1)
  WAIC2    = -2*(lppd - q2)
  mat=matrix(nrow=1, ncol=3)
  mat[1,]=c(WAIC1, WAIC2, DIC)
  colnames(mat) = c('WAIC1', 'WAIC2', 'DIC')
  mat
}