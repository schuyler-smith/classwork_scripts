calc_xi2 <- function(E, K, l_val){
  xi_2 <- rep(0, l_val)
  for(l in 0:(l_val-1)){
    xi_2[l+1] <- xi_2[l+1] + sum(E[,seq(l+1,l_val*(K-1),l_val)])
  }
  xi_2 <- xi_2/sum(xi_2)
  return(xi_2)
}