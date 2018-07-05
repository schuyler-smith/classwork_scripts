
simulate_MC_reads <- function(
  MSA, 
  fasta_file = MSA@fasta_file , 
  order = MSA@order ,
  structures = MSA@structures ,
  alpha = MSA@alpha , 
  P.hat = MSA@transition_frequencies ,
  likelihoods=MSA@likelihoods ){
  
number_structures <- (width(readBStringSet(fasta_file)))

priors <- unique(P.hat[[1]])
P.hat <- data.matrix(dcast(P.hat, prior_state ~ transition_state, value.var = "probability", fill = 0))[,-1]
rownames(P.hat) <- priors

structures <- colnames(P.hat)
simulated_sequences <- character(length=length(number_structures))

print(paste("simulating reads for order ",order,sep=""))
for(x in 1:length(number_structures)){
  progress(x, length(number_structures))
  Sys.sleep(0.001)
  
  if(number_structures[x] < order){
    start <- sample(rownames(alpha), size = 1, prob = alpha$alpha)
    start <- strtrim(start,number_structures[x])
    initial <- paste(start,"*",sep="")
    n <- sum(alpha[grepl(glob2rx(initial), rownames(alpha)),])
    simulated_sequences[x] <- start
  }

  else if(number_structures[x] == order){
    start <- sample(rownames(alpha), size = 1, prob = alpha$alpha)
    simulated_sequences[x] <- start
  }

  else{
    chain <- character((number_structures[x]-(order-1)))
      while(chain[length(chain)] == ""){
        start <- sample(rownames(alpha), size = 1, prob = alpha$alpha)
        chain[1] <- start
        try(for(y in 1:(number_structures[x]-order)){
          next_in_chain <- sample(structures, size = 1, prob = P.hat[start,])
          chain[(y+1)] <- next_in_chain
          start <- paste(substr(start,2,order),next_in_chain,sep="")
          }, silent = TRUE)}
      simulated_sequences[x] <- paste(chain, collapse = "")
    }
  }
return(simulated_sequences)
}
