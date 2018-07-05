calc_gamma <- function(M, structures, reference, biallelic_regions, gamma, E, K, l_val, IP){
  which_gamma<-function(x,y){return(((x-1)*(4))+y)}
  gamma <- matrix(0, length(M), length(structures)^2)
  for(i in 1:length(M)){
    m <- unlist(strsplit(M[i],""))
    for(k in 0:(K-2)){
      for(l in 0:(l_val-1)){
        for(j in IP){
          if(m[j+k] != reference[j]){
            id <- which_gamma(which(structures == reference[j]),which(structures == m[j+k]))
            gamma[i, id] <- gamma[i, id]+1}}
        for(j in biallelic_regions){
          if(reference[j] == "Y"){
            if(l < 2 && m[j+k] != "C"){
                id <- which_gamma(which(structures == "C"),which(structures == m[j+k]))
                gamma[i, id] <- gamma[i, id]+1}
            if(l >= 2 &&  m[j+k] != "T"){
                id <- which_gamma(which(structures == "T"),which(structures == m[j+k]))
                gamma[i, id] <- gamma[i, id]+1}
            }
          else if(reference[j] == "R"){  
            if(l %% 2 == 0 && m[j+k] != "A"){
              id <- which_gamma(which(structures == "A"),which(structures == m[j+k]))
              gamma[i, id] <- gamma[i, id]+1}
            if(l %% 2 != 0 && m[j+k] != "G"){
              id <- which_gamma(which(structures == "G"),which(structures == m[j+k]))
              gamma[i, id] <- gamma[i, id]+1}
          }
        }
        gamma[i,] <- gamma[i,]*E[i, ((k*(K-1))+l)+1]
      }
    }
  }
  gamma <- apply(gamma,2,sum)
  for(l in seq(1, length(gamma), l_val)){
    gamma[(l:(l+l_val-1))] <- gamma[(l:(l+l_val-1))]/sum(gamma[(l:(l+l_val-1))])
  }
  return(gamma)
}