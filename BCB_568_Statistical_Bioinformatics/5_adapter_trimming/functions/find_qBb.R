calc_qBb <- function(M, K, l_val, structures, E, N_regions){
  nucleotide_E <- matrix(0, length(M), length(structures))
  for(i in 1:length(M)){
    m <- unlist(strsplit(M[i],""))
    for(k in 0:(K-2)){
      n_regions <- N_regions+k
      nucleotide_counts <- apply(array(structures),1,FUN=function(b){
        sum(m[0:k] %in% b, m[n_regions] %in% b)})
      for(l in 0:(l_val-1)){
        nucleotide_E[i,] <- nucleotide_E[i,] + (nucleotide_counts * E[i, ((k*(K-1))+l)+1])
      }
    }
  }
  
  q_Bb <- apply(nucleotide_E,2,sum)/sum(apply(nucleotide_E,2,sum))
return(q_Bb)}  

