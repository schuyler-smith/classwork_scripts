
calc_delta <- function(M, reference, K, l_val, IP, E, biallelic_regions){
  delta_A <- (rep(0, length(reference)))
  delta_B <- (rep(0, length(reference)))
  
  for(i in 1:length(M)){
    m <- unlist(strsplit(M[i],""))
    for(k in 0:(K-2)){
      for(l in 0:(l_val-1)){
        delta_Aikl <- (rep(0, length(reference)))
        delta_Bikl <- (rep(0, length(reference)))
        for(j in (IP)){
          if(m[j+k] == reference[j]){
            delta_Aikl[j] <- 1}
          else{
            delta_Bikl[j] <- 1}}
        for(j in biallelic_regions){
          if(reference[j] == "Y"){
            if(l < 2){
              if(m[j+k] == "C"){
                delta_Aikl[j+k] <- 1}}
            if(l >= 2){
              if(m[j+k] == "T"){
                delta_Aikl[j+k] <- 1}}
            else{delta_Bikl[j+k] <- 1}}
          else if(reference[j] == "R"){  
            if(l %% 2 == 0){
              if(m[j+k] == "A"){
                delta_Aikl[j+k] <- 1}}
            if(l %% 2 != 0){
              if(m[j+k] == "G"){
                delta_Aikl[j+k] <- 1}}
            else{delta_Bikl[j+k] <- 1}}
        }
      delta_A <- delta_A + delta_Aikl*E[i, ((k*(K-1))+l)+1]
      delta_B <- delta_B + delta_Bikl*E[i, ((k*(K-1))+l)+1]
      }
    }
  }
  return(delta_A/(delta_A+delta_B))
}



