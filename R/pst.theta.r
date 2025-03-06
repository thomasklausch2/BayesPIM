pst.theta = function(C, a){
  x = table(C)
  alpha = a + x
  rdirichlet(1, alpha)
}
