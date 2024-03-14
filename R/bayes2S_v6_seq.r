#' bayes.2S_seq
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
bayes.2S_seq = function(Vobs, kappa, Z.X = NULL, Z.W = NULL, r = rep(0, length(Vobs)), parallel = T, update.kappa = F, prev = T,
                        ndraws=1000, chains=3, thining=1, prop.sd.X=NULL, naive.run.prop.sd.X = prop.sd.X,
                        beta.prior.X = 4, sig.prior.X = 10, tau.w = 1, k.prior = sqrt(10), fix.sigma.X = F,
                        prev.run = NULL, dist.X = 'weibull', update.burnin = T,
                        beta.prior = 't', vanilla = F, ndraws.naive = 5e3, update.till.converge = F, ndraws.update = NULL, maxit=Inf,
                        conv.crit = 'upper', kappa.prior = NULL, fix.k = F, par.exp=T, collapsed.g = F){
  t0 = Sys.time()
  if(!is.numeric(maxit)) stop('maxit has to be numeric.')
  
  if(!is.null(prev.run)){
    chains = length(prev.run$par.X.all)
    dims = dim(as.matrix(prev.run$par.X.all[1]))
    start.val.X = matrix(ncol=dims[2], nrow= chains )
    for(i in 1:chains){
      start.val.X[i,] = as.matrix(prev.run$par.X.all[i])[dims[1],]
    }
    Z.X = prev.run$Z.X
    Z.W = prev.run$Z.W
    r   = prev.run$r
    Vobs = prev.run$Vobs
    kappa = prev.run$kappa[length(prev.run$kappa)]
    update.kappa = prev.run$update.kappa
    prev = prev.run$prev
    kappa.prior = prev.run$kappa.prior
    if(is.infinite(maxit)) maxit = prev.run$maxit
    vanilla = prev.run$vanilla
    par.exp = prev.run$par.exp
    ndraws.prev = prev.run$ndraws
    if(is.null(ndraws.update)) ndraws = prev.run$ndraws else ndraws = ndraws.update
    if(is.null(prop.sd.X)){ prop.sd.X =  prev.run$prop.sd.X }
    beta.prior.X = prev.run$priors$beta.prior.X
    sig.prior.X = prev.run$priors$sig.prior.X
    thining = prev.run$thining
    X.prev = prev.run$X
    C.prev = prev.run$C
    dist.X = prev.run$dist.X
    fix.sigma.X = prev.run$fix.sigma.X
    beta.prior = prev.run$beta.prior
    kappa.prior = prev.run$kappa.prior
    fix.k = prev.run$fix.k
    prev.runtime = prev.run$runtime
    cat('Updating previous MCMC run. \n')
  }
  
  if(update.kappa & !is.null(kappa.prior)){
    f = find.ab(m=kappa.prior[1], s=kappa.prior[2])
    kappa.ab = c(f$a,f$b)
    if(f$conv == 1 | f$value>0.01) stop('Did not find hyperparameters for pi(kappa|a,b). Specify different kappa.prior.')
  }
  if(update.kappa & is.null(kappa.prior)) kappa.ab = c(1,1)
  
  # Prepare data
  g.fixed = sapply(Vobs, function(x) (length(x) == 1))
  for(i in 1:length(Vobs)) if(g.fixed[i]) Vobs[[i]] = c(0, Inf)
  
  burnin=round(ndraws/2)
  pobs = P_vobs (Vobs, kappa = kappa)
  
  d = sapply(Vobs, function(x) as.numeric(is.finite(x[length(x)]) )+1)
  L = sapply(Vobs, function(x) x[length(x)-1])
  R = sapply(Vobs, function(x) x[length(x)])
  
  pobs_vec = unlist(pobs)
  Vobs_L = unlist( lapply(Vobs, function(x) x[1:(length(x)-1)] ) )
  Vobs_R = unlist( lapply(Vobs, function(x) x[2:(length(x))] ) )
  m   = sapply(Vobs, length)-1
  
  #
  
  if(is.null(d) | is.null(L) |is.null(R)){
    stop("Provide data.")
  }
  if(is.null(prop.sd.X)  ){
    stop("Provide proposal SD")
  }
  
  if(!is.null(Z.X)) {
    Z.X = as.matrix(Z.X)
    if( is.null(colnames(Z.X)) ) colnames(Z.X) = paste('ZX.',1:ncol(Z.X))
    p.X = ncol(Z.X)
    Z1.X = cbind(1,Z.X)
  }
  if(is.null(Z.X)){
    p.X = 0
    Z1.X  = matrix(1, ncol=1, nrow=length(d))
  }
  
  if(!is.null(Z.W)) {
    Z.W = as.matrix(Z.W)
    if( is.null(colnames(Z.W)) ) colnames(Z.W) = paste('ZW.',1:ncol(Z.W))
    Z1.W = cbind(1,Z.W)
  }
  if(is.null(Z.W)){
    Z1.W  = matrix(1, ncol=1, nrow=length(d))
  }
  
  p1.X = ncol(Z1.X)
  p1.W = ncol(Z1.W)
  
  sig_inv = solve(t(Z1.W) %*% Z1.W + diag(rep(tau.w^-1,ncol(Z1.W)) ))
  sig_inv_Xt = sig_inv %*% t(Z1.W)
  
  # Get random seeds
  s = sample(1:10^5,chains, replace=F)
  
  if(!vanilla & is.null(prev.run)){
    adapted = 2
    cat('Searching starting values by one naive run \n')
    done = F
    while(!done){
      mod.simple = bayes.2S(Vobs[!g.fixed], kappa, Z.X = Z.X[!g.fixed,], r=r[!g.fixed], parallel = T,
                            ndraws=ndraws.naive, chains=3, thining=1, prop.sd.X=naive.run.prop.sd.X,
                            beta.prior.X = 4, sig.prior.X = 10, fix.sigma.X = fix.sigma.X,
                            prev.run = NULL, dist.X = dist.X, update.burnin = T,
                            beta.prior = 't', vanilla = T)
      pars.simple = as.matrix(mod.simple$par.X.bi)
      ini = pars.simple[sample(1: nrow(pars.simple), chains, replace=F),]
      ac.rates = apply(mod.simple$ac.X,2, mean)
      if( sum(ac.rates > .85) >0) {
        cat('Naive run acceptance rates >0.85. Increasing naive.run.proposal.sd.X.\n')
        naive.run.prop.sd.X = naive.run.prop.sd.X*adapted
        adapted = adapted + 1
      }
      if( sum(ac.rates < .15) >0) {
        cat('Naive run acceptance rates <0.15. Decreasing naive.run.proposal.sd.X.\n')
        naive.run.prop.sd.X = naive.run.prop.sd.X/adapted
        adapted = adapted + 1
      }
      if(adapted>20) stop('Cannot find naive.run.prop.sd.X that produces good acceptance rate.')
      if(sum(ac.rates > .85) == 0 & sum(ac.rates < .15) == 0) done = T
    }
    cat('Now doing main run. \n')
  }
  
  ## Start run
  cat(paste("Starting Gibbs sampler with", chains, "chains and", ndraws, "iterations.\n"))
  run = out = list()
  for(j in 1:chains){
                
                  set.seed(s[j])
                  
                  # Begin Gibbs
                  ac.X = ac.S = 1
                  
                  # Init
                  if(!vanilla & is.null(prev.run)){
                    beta.X.ini  = ini[j, 2:(1+p.X)]
                    mu.X.ini    = ini[j, 1]
                    sigma.X.ini = ini[j, ncol(ini)]
                    beta_w.ini  = runif(p1.W,-1,1)
                    beta_w      = matrix(NA, nrow = ndraws +1, ncol = p1.W)
                    beta_w[1,]  = beta_w.ini
                    C.aug       = rep(0, length(L))
                    C.aug[g.fixed] = 1
                    if(update.kappa) kappa  = runif(1, 0.2,0.8)
                    if(dist.X == 'gengamma') log.k.ini = runif(1, -2, 2)
                  } else{
                    beta.X.ini  = runif(p.X, -1,1)
                    mu.X.ini    = runif(1, -1,1)
                    sigma.X.ini = runif(1, -2,2)
                    if(dist.X == 'gengamma') log.k.ini = runif(1, -2, 2)
                  }
                  
                  if(fix.sigma.X) sigma.X.ini = log(sig.prior.X)
                  if(fix.k) log.k.ini = log(k.prior)
                  
                  if(!is.null(prev.run)){
                    mu.X.ini = start.val.X[j,1]
                    if(p.X > 0) beta.X.ini = start.val.X[j,2:(p.X+1)] else beta.X.ini = NULL
                    sigma.X.ini = log(start.val.X[j,(p.X+2)])
                    
                    if(!vanilla){
                      beta_w      = matrix(NA, nrow = ndraws +1, ncol = p1.W)
                      beta_w.ini = start.val.X[j,(p.X+3):(p.X+p1.W+2)]
                      beta_w[1,]  = beta_w.ini
                    }
                    
                    n.prev = length(Vobs)
                    X = X.prev[ ((n.prev*(j-1))+1): (n.prev*j) ]
                    if(!vanilla){
                      C.aug = C.prev[ ((n.prev*(j-1))+1): (n.prev*j) ]
                    }
                    if(update.kappa) kappa = start.val.X[j,ncol(start.val.X)]
                  }
                  
                  if(dist.X != 'gengamma'){
                    cur.par.Xreg = matrix(ncol=p.X+2, nrow=ndraws+1)
                    cur.par.Xreg[1,1:(p.X+1)] = c(mu.X.ini,beta.X.ini)
                    cur.par.Xreg[1,(p.X+2)]   = sigma.X.ini
                  }    else{
                    cur.par.Xreg = matrix(ncol=p.X+3, nrow=ndraws+1)
                    cur.par.Xreg[1,1:(p.X+1)] = c(mu.X.ini,beta.X.ini)
                    cur.par.Xreg[1,(p.X+2)]   = sigma.X.ini
                    cur.par.Xreg[1,(p.X+3)]   = log.k.ini
                  }
                  
                  if(is.null(prev.run)){
                    X   = pst.X.2S( cbind(rep(.1,length(d)),NA) , rgamma(length(d),.1), d = d, L= L, R= R, dist = "exp")
                    if(dist.X == 'weibull' | dist.X == 'loglog') cur.par.X = trans.par(Z1.X, par = cur.par.Xreg[1,])
                    if(dist.X == 'lognormal') cur.par.X = trans.par.ind.norm(Z1 = Z1.X, p = cur.par.Xreg[1,1:p1.X], v= cur.par.Xreg[1,(p1.X+1)])
                    if(dist.X == 'gengamma')  cur.par.X = trans.par.gengamma(Z1 = Z1.X, par = cur.par.Xreg[1,])
                    X   = pst.X.2S(par=cur.par.X, d = d, L= L, R= R, dist = dist.X)
                  }
                  
                  i=1
                  log.pst.X  = pst.aft(par=cur.par.Xreg[i,, drop=F], t=X, Z=Z1.X, tau= beta.prior.X, sig.prior=sig.prior.X,
                                       k.prior = k.prior,
                                       dist=dist.X, beta.prior = beta.prior)
                  if(is.infinite(log.pst.X) ) stop('Bad starting values')
                  
                  
                  prop.sd.X.mat = diag(prop.sd.X^2,p1.X+1)
                  if(dist.X == 'gengamma') prop.sd.X.mat = diag(prop.sd.X^2,p1.X+2)
                  
                  for(i in 1:ndraws){
                    # if (i %% 100 == 0) {
                    #   cat("Completed", i, "of", n, "iterations\n")
                    # }
                    #Update  X parameters
                    # if(i == burnin) {
                    #   S = cov(cur.par.Xreg[(burnin-101):(burnin-1),])
                    #   prop.sd.X.mat = 2.4^2/(nrow(Z1.X)+p1.X+1) * S #+ 0.0001 * diag(1, p1.X+1)
                    # }
                    mh  = mhstep.aft(x=cur.par.Xreg[i,, drop=F], t=X, Z=Z1.X, tau= beta.prior.X, sig.prior=sig.prior.X, k.prior = k.prior,
                                     prop.var=prop.sd.X.mat, dist=dist.X, fix.sigma=fix.sigma.X, fix.k = fix.k)
                    cur.par.Xreg[i+1,] = mh$s
                    
                    # Assess acceptance for MH steps
                    ac.X[i+1] = cur.par.Xreg[i+1,1] != cur.par.Xreg[i,1]
                    
                    # augment X
                    if(dist.X == 'weibull' | dist.X == 'loglog') cur.par.X = trans.par(Z1.X, par = cur.par.Xreg[(i+1),])
                    if(dist.X == 'lognormal') cur.par.X = trans.par.ind.norm(Z1 = Z1.X, p = cur.par.Xreg[(i+1),1:p1.X], v= cur.par.Xreg[(i+1),(p1.X+1)])
                    if(dist.X == 'gengamma')  cur.par.X = trans.par.gengamma(Z1 = Z1.X, par = cur.par.Xreg[(i+1),])
                    
                    if(!vanilla & prev) {
                      # Update pobs
                      if(update.kappa) pobs = P_vobs_Rcpp(Vobs, kappa = kappa[i])
                      
                      # Update X
                      if(dist.X %in% c('weibull','lognormal')) aug.X = augment.X_rcpp( unlist(pobs), Vobs, Vobs_L, Vobs_R, cur.par.X, dist.X, C=C.aug, collapsed.g = collapsed.g)
                      else aug.X = augment.X( unlist(pobs), Vobs, Vobs_L, Vobs_R, cur.par.X, dist.X, C=C.aug, collapsed.g = collapsed.g)
                      X = aug.X$X
                      
                      # Calculate theta1
                      mu_w    = as.numeric(Z1.W %*% as.matrix(as.numeric(beta_w[i,])))
                      theta1  = pnorm( mu_w )
                      
                      # Update C.aug
                      if(!collapsed.g) C.aug   = augment_C_Rcpp(pobs, Vobs, X, kappa = kappa[i], theta1 = theta1,
                                                                r= r, g_fixed = g.fixed)
                      else  C.aug   = augment_C_collapsed_rcpp(w_sums = aug.X$sums, Vobs, kappa = kappa[i], theta1 = theta1, r= r, g_fixed = g.fixed)
                      
                      W.aug   = augment.W(g = C.aug, mu_w )
                      if(par.exp){
                        alpha_sq     = fc_w_par.exp_Haar(y = W.aug, X = Z1.W, sig_inv_Xt = sig_inv_Xt) #solve(t(X)X) and t(X) can be calculated outside loop
                        beta_w[i+1,] = fc_beta(X = Z1.W, W.aug/sqrt(alpha_sq), sig_inv_Xt = sig_inv_Xt, sig_inv=sig_inv)
                      } else{
                        beta_w[i+1,] = fc_beta(X = Z1.W, W.aug, sig_inv_Xt = sig_inv_Xt, sig_inv=sig_inv)
                      }
                      
                      if(update.kappa) kappa[i+1] = fc_kappa_rcpp(Vobs, j_ = aug.X$k, a=kappa.ab[1], b=kappa.ab[2], g=C.aug, r = r, g_fixed= g.fixed)
                      else kappa[i+1] = kappa[i]
                    }
                    if(!vanilla & !prev) {
                      if(update.kappa) pobs = P_vobs (Vobs, kappa = kappa[i]) #
                      aug.X = augment.X( unlist(pobs), Vobs, Vobs_L, Vobs_R, cur.par.X, dist.X, C=C.aug)
                      X = aug.X$X
                      if(update.kappa) kappa[i+1] = pst.kappa.noprev( Vobs, j_ = aug.X$k, a=kappa.ab[1], b=kappa.ab[2])
                      else kappa[i+1] = kappa[i]
                    }
                    if(vanilla)  X = pst.X.2S(par=cur.par.X, d = d, L= L, R= R, dist = dist.X)
                    X[X==0] = 10^-300
                  }
                  out=list()
                  out$X = X
                  if(!vanilla & prev) {
                    out$par.X = cbind(cur.par.Xreg[-1,], beta_w[-1,])
                    out$C.aug = C.aug
                  }
                  if(!vanilla & !prev) {
                    out$par.X = cur.par.Xreg[-1,]
                    out$C.aug = C.aug
                  }
                  if(vanilla)  out$par.X = cur.par.Xreg[-1,]
                  out$kappa = kappa[-1]
                  out$ac.X = ac.X
                  run[[j]] = out
  } 
  
  # Unwrap chains wihout trimming
  mcmc.par.X= list()
  i=1
  par.X = run[[i]]$par.X
  if(vanilla) {
    if(dist.X != 'gengamma'){
      colnames(par.X) = c("Intercept", colnames(Z.X), "sigma")
      par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
    }
    if(dist.X == 'gengamma'){
      colnames(par.X) = c("Intercept", colnames(Z.X), "sigma",'k')
      par.X[,ncol(par.X)-1] = exp(par.X[,ncol(par.X)-1])
      par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
    }
  }
  if(!vanilla & prev) {
    colnames(par.X) = c("Intercept X", colnames(Z.X), "sigma", "Intercept W", colnames(Z.W))
    par.X[,ncol(Z.X)+2] = exp(par.X[,ncol(Z.X)+2])
  }
  if(!vanilla & !prev) {
    colnames(par.X) = c("Intercept", colnames(Z.X), "sigma")
    par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
  }
  if(update.kappa) {
    par.X = cbind(par.X, run[[i]]$kappa)
    colnames(par.X)[ncol(par.X)] = 'kappa'
  }
  mcmc.par.X[[i]] = mcmc(par.X)
  ac.X = run[[i]]$ac.X
  X.draw = run[[i]]$X
  if(!vanilla) C.draw = run[[i]]$C.aug
  
  if(length(run)>1){
    for(i in 2:length(run)){
      par.X = run[[i]]$par.X
      if(vanilla) {
        if(dist.X != 'gengamma'){
          colnames(par.X) = c("Intercept", colnames(Z.X), "sigma")
          par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
        }
        if(dist.X == 'gengamma'){
          colnames(par.X) = c("Intercept", colnames(Z.X), "sigma",'k')
          par.X[,ncol(par.X)-1] = exp(par.X[,ncol(par.X)-1])
          par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
        }
      }
      if(!vanilla & prev) {
        colnames(par.X) = c("Intercept X", colnames(Z.X), "sigma", "Intercept W", colnames(Z.W))
        par.X[,ncol(Z.X)+2] = exp(par.X[,ncol(Z.X)+2])
      }
      if(!vanilla & !prev) {
        colnames(par.X) = c("Intercept", colnames(Z.X), "sigma")
        par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
      }
      if(update.kappa) {
        par.X = cbind(par.X, run[[i]]$kappa)
        colnames(par.X)[ncol(par.X)] = 'kappa'
      }
      mcmc.par.X[[i]] = mcmc(par.X)
      ac.X = cbind(ac.X, run[[i]]$ac.X)
      X.draw = c(X.draw, run[[i]]$X)
      if(!vanilla) C.draw = c(C.draw, run[[i]]$C.aug)
    }}
  
  par.X.all = mcmc.list(mcmc.par.X)
  
  if(!is.null(prev.run)){
    par.X.all = bind.mcmclists(prev.run$par.X.all, par.X.all)
    ac.X.cur = ac.X
    ac.X = rbind(prev.run$ac.X, ac.X)
    mc.prev  = nrow(as.matrix(prev.run$par.X.all[1]))
    if(update.burnin) burnin   = round((mc.prev + ndraws)/2)
  }
  nr       = nrow(as.matrix(par.X.all[1]))
  par.X.bi = trim.mcmc(par.X.all, burnin = burnin, thining = thining)
  #
  # Undo Vobs recoding
  for(i in 1:length(Vobs)) if(Vobs[[i]][1] ==0 & is.infinite(Vobs[[i]][2])) Vobs[[i]] = 0
  
  
  t1 = Sys.time()
  runtime = t1-t0
  if(!is.null(prev.run)) runtime = runtime + prev.runtime
  # ## save built info
  dat = data.frame(L=L, R=R)
  priors = list()
  priors$beta.prior.X = beta.prior.X
  priors$sig.prior.X = sig.prior.X
  
  
  ret = list()
  ret$par.X.all = par.X.all
  ret$par.X.bi = par.X.bi
  ret$X = X.draw
  ret$ac.X = ac.X
  if(!is.null(prev.run)){
    ret$ac.X.cur = ac.X.cur
  }
  ret$dat = dat
  ret$Vobs = Vobs
  ret$Z.X = Z.X
  ret$Z.W = Z.W
  ret$r    = r
  ret$kappa = kappa
  ret$kappa.prior = kappa.prior
  ret$maxit = maxit
  ret$prev = prev
  ret$vanilla = vanilla
  ret$par.exp = par.exp
  ret$update.kappa = update.kappa
  ret$priors = priors
  ret$thining = thining
  ret$prop.sd.X = prop.sd.X
  ret$dist.X = dist.X
  ret$fix.sigma.X = fix.sigma.X
  ret$burnin = burnin
  ret$beta.prior = beta.prior
  ret$runtime = runtime
  ret$ndraws = ndraws
  if(!vanilla) ret$C = C.draw
  ret$kappa.prior = kappa.prior
  ret$fix.k = fix.k
  
  if(update.till.converge){
    if(vanilla) {
      g <- gelman.diag(ret$par.X.bi[,-(p1.X+2)])
      effs <- effectiveSize( ret$par.X.bi[,-(p1.X+2)] )
    }
    else{
      g <- gelman.diag(ret$par.X.bi)
      effs <- effectiveSize( ret$par.X.bi )
    }
    
    if(conv.crit == 'point') ind = 1
    if(conv.crit == 'upper') ind = 2
    cri1 = sum(g[[1]][,ind]<1.1) == length(g[[1]][,ind])
    cri2 = sum(effs > chains * 2 * 5) == length(effs)
    cri3 = dim(ret$par.X.all[[1]])[1] >= maxit
    
    while(!(cri1 & cri2) & !cri3){
      cat(paste('Completed', dim(ret$par.X.all[[1]])[1], 'draws. \n'))
      cat(paste('Acceptance rates were:', paste(round( apply(ret$ac.X,2,mean), 2), collapse = ', '),'\n'))
      cat('Not converged. \n')
      # cat(paste( round(g[[1]][,ind], 2), '\n'))
      # cat(paste( round(effs), '\n'))
      if(is.null(ndraws.update)) ret = bayes.2S(prev.run = ret, ndraws = ndraws)
      else ret = bayes.2S(prev.run = ret, ndraws.update = ndraws.update)
      if(vanilla) {
        g <- gelman.diag(ret$par.X.bi[,-(p1.X+2)])
        effs <- effectiveSize( ret$par.X.bi[,-(p1.X+2)] )
      }
      else{
        g <- gelman.diag(ret$par.X.bi)
        effs <- effectiveSize( ret$par.X.bi )
      }
      cri1 = sum(g[[1]][,ind]<1.1) == length(g[[1]][,ind])
      cri2 = sum(effs > chains * 2 * 5) == length(effs)
      cri3 = dim(ret$par.X.all[[1]])[1] >= maxit
    }
    cat(paste('Completed', dim(ret$par.X.all[[1]])[1], 'draws. \n'))
    if(!cri3) cat('Converged. \n') else cat('Maxit reached. Not converged. \n')
    ret$convergence = !cri3
  }
  return(ret)
}
