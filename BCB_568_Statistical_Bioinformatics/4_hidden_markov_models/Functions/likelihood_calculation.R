
likelihood_calculation <- function(MSA, 
                                   fasta_file = MSA@fasta_file, 
                                   order = MSA@order,
                                   structures = MSA@structures,
                                   alpha = MSA@alpha, 
                                   transition_frequencies = MSA@transition_frequencies){
  
  sequences <- paste(readBStringSet(fasta_file))
  number_structures <- width(sequences)
  likelihood <- numeric(length=length(number_structures))
  initials <- apply(array(sequences), 1, FUN=function(y){strtrim(y,order)})
  
  print("calculating likelihoods for reads")
  for(x in 1:length(number_structures)){
          progress(x, length(number_structures))
          if(x == length(number_structures)) cat("Done!\n")
    
    transitions <- unlist(strsplit(sequences[x],""))[-c(1:order)]
    start <- initials[x]
    
    if(number_structures[x] < order){
      n <- 0
      while(n == 0){
        initial <- paste(start,"*",sep="")
        start <- strtrim(start,(nchar(start)-1))
        n <- sum(alpha[grepl(glob2rx(initial), rownames(alpha)),])}
      likelihood[x] <- log(n*((1/length(structures))^(order-(nchar(initial)-1))))
    }
    
    else if(number_structures[x] == order){
      likelihood[x] <- log(alpha[initials[x],])
    }
    
    else{
      likelihood[x] <- sum(log(alpha[initials[x],]), 
        apply(array(1:(number_structures[x]-order)),1,FUN=function(y){
          next_in_chain <- transitions[y]
          transition_likelihood <- log(transition_frequencies[start,next_in_chain])
          start <<- paste(substr(start,2,order),next_in_chain,sep="")
          return(transition_likelihood)}))
    }
  }
  return(likelihood)
}
