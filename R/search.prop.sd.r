#' search.prop.sd
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
search.prop.sd = function(m, ndraws = 1000, succ.min = 3, acc.bounds.X =c(0.2,0.25)){
  found.X = found.S = F
  it = 1; succ = 0;
  while(succ!=succ.min){
    cat(paste('Iteration',it,'\n'))
    if(it == 1) { ac.X.cur = mean(m$ac.X)
    prop.sd.X = m$prop.sd.X}
    if(it >1) {
      sink(file = "NUL", type = "output")
      m = bayes.2S( prev.run = m, ndraws.update = ndraws, prop.sd.X = prop.sd.X)
      sink(type = "output")
      ac.X.cur = mean(m$ac.X.cur)
    }
    acc.bounds.mean = (acc.bounds.X[2]-acc.bounds.X[1])/2 + acc.bounds.X[1]
    ac.X = mean(m$ac.X)
    cat( paste('Acceptance rate was:', round(ac.X.cur, 3),'\n' ))
    ac.X = ac.X.cur
    if((ac.X > acc.bounds.X[1] & ac.X < acc.bounds.X[2]) ){
      found.X = T} else{
      if(ac.X < acc.bounds.X[1]) {
        dif =  1-(acc.bounds.X[1] - ac.X)/acc.bounds.X[1]
        prop.sd.X=prop.sd.X * dif
      }
      if(ac.X > acc.bounds.X[2]) {
        dif =  1+(ac.X - acc.bounds.X[2])/acc.bounds.X[2]
        prop.sd.X=prop.sd.X * dif
      }
      found.X = F
      cat(paste('prop.sd.X is set to', round(prop.sd.X,3),'\n'))
    }
    it=it+1
    if(found.X ){
      ndraws = ndraws*2
      succ = succ +1
      found.X = F
      if(succ!=succ.min) cat(paste('Success. Doubling number of MCMC draws:',ndraws,'\n'))
      }
  }
  cat('Finished calibrating proposal variance. \n')
  ret= list()
  ret$prop.sd.X = prop.sd.X
  ret$ac.X      = ac.X
  #ret$mod       = m
  ret
}
