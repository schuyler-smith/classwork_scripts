calc_E <- function(M, 
                  reference = reference, 
                  K = K, 
                  l_val = l_val, 
                  q_Sb = q_Sb, 
                  q_Bb = q_Bb, 
                  gamma = q_Bb, 
                  delta = delta, 
                  xi_1 = xi_1, 
                  xi_2 = xi_2,
                  structures = structures,
                  IP = IP
                 ){
  which_gamma<-function(x,y){return(((x-1)*(4))+y)}
  biallelic_regions = which((unlist(strsplit(reference,"")) %in% c("Y","R")))
  N_regions = which((unlist(strsplit(reference,""))) == "N")
  log_likelihoods <- matrix(0,length(M),((K-1)*l_val)+1)
  
  for(i in 1:length(M)){
    m <- unlist(strsplit(M[i],""))
    for(k in 0:(K-2)){
      log_xi_k <- log(xi_1[k+1])
      sequence_regions <- (k+length(reference)+1):length(m)
      for(l in 0:(l_val-1)){
        log_xi_l <- log(xi_2[l+1])
        some_sums <- sum(log(q_Bb[apply(array(unlist(strsplit(M[i],""))[1:(k+1)]), 1, FUN=function(x){which(structures == x)})]))
        for(j in IP){
          if(m[j+k] == reference[j]){some_sums <- some_sums + log(delta[j])}
          else{some_sums <- some_sums + log(1-delta[j]) + log(gamma[which_gamma(which(structures == reference[j]), which(structures == m[j+k]))])}}
        for(j in c(N_regions,sequence_regions)){
          some_sums <- some_sums + log(q_Bb[structures == m[j]])}
        
        for(j in biallelic_regions){
          if(reference[j] == "Y"){
            if(l < 2){
              if(m[j+k] == "C"){some_sums <- some_sums + log(delta[j+k])}
              else{some_sums <- some_sums + (log(1-delta[j+k])+ log(gamma[which_gamma(which(structures =="C"),which(structures == m[j+k]))]))}}
            if(l >= 2){
              if(m[j+k] == "T"){some_sums <- some_sums + log(delta[j+k])}
              else{some_sums <- some_sums + (log(1-delta[j+k])+ log(gamma[which_gamma(which(structures == "T"),which(structures == m[j+k]))]))}}
          }
          else if(reference[j] == "R"){  
            if(l %% 2 == 0){
              if(m[j+k] == "A"){some_sums <- some_sums + log(delta[j+k])}
              else{some_sums <- some_sums + (log(1-delta[j+k])+ log(gamma[which_gamma(which(structures =="A"),which(structures == m[j+k]))]))}}
            if(l %% 2 != 0){
              if(m[j+k] == "G"){some_sums <- some_sums + log(delta[j+k])}
              else{some_sums <- some_sums + (log(1-delta[j+k])+ log(gamma[which_gamma(which(structures =="G"),which(structures == m[j+k]))]))}}
          }
        }
        log_likelihoods[i, ((k*(K-1))+l)+1] <- some_sums + log_xi_k + log_xi_l
      }
    }
  }

  for(i in 1:length(M)){
    m <- unlist(strsplit(M[i],""))
    N <- apply(array(structures),1,FUN=function(b){sum(m %in% b)})
    N_table <- data.table(n = N, q_Sb = log(q_Sb))
    log_likelihoods[i,((K-1)*l_val)+1] <- sum(apply(N_table, 1, prod), xi_1[K])
  }
  
  E <- t(apply(log_likelihoods, 1, FUN=function(x){
    ll_max <- max(x)
    probibilities <- apply(array(x),1,FUN=function(y){
      10^(y-ll_max)})
    probibilities/sum(probibilities)
  }))

return(E)}






