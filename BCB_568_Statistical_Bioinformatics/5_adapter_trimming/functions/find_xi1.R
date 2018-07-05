calc_xi1 <- function(M, E, K, l_val){
  xi_1 <- rep(0, K)
  for(k in 0:(K-2)){
    xi_1[k+1] <- xi_1[k+1] + sum(E[,(((k*(K-1)))+1):((((k*(K-1)))+1)+l_val-1)])
  }
  xi_1[K] <- sum(E[,((K-1)*l_val)+1])
  xi_1 <- xi_1/length(M)
  return(xi_1)
}