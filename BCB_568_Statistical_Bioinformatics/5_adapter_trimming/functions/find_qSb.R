calc_qSb <- function(M, K, l_val, structures, E, reference){
  nucleotide_E <- matrix(0, length(M), length(structures))
  for(i in 1:length(M)){
    m <- unlist(strsplit(M[i],""))
    for(k in 0:(K-2)){
      mk <- m[(k+length(reference)+1):length(m)]
      nucleotide_counts <- apply(array(structures),1,FUN=function(b){
        sum(mk %in% b)})
      for(l in 0:(l_val-1)){
        nucleotide_E[i,] <- nucleotide_E[i,] + (nucleotide_counts * E[i, ((k*(K-1))+l)+1])}
      if(k == 0){
        nucleotide_counts <- apply(array(structures),1,FUN=function(b){
          sum(m %in% b)})
        nucleotide_E[i, ] <- nucleotide_E[i,] + (nucleotide_counts*E[i, ((K-1)*l_val)+1])}
    }
  }
  q_Sb <- apply(nucleotide_E,2,sum)/sum(apply(nucleotide_E,2,sum))
  return(q_Sb)
}
