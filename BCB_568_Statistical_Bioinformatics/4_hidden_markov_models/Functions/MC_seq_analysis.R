

MC_seq_analysis <- function(fasta_file = fasta_file, order = 1){
                      setClass("MC_Sequence_Analysis", representation(
                                                fasta_file="character",
                                                order="numeric",
                                                structures="factor", 
                                                alpha="data.frame",
                                                transition_frequencies="data.table",
                                                likelihood="numeric"
                                                )
         )
transitions <- data.table(zwc::fasta_wc(k=order+1, fasta_file))
seq_data <- readBStringSet(fasta_file)
seqs <-paste(seq_data)  

transitions <- transitions[order(word)]
structures <- sort(unique(unlist(apply(array(transitions[[1]]),1,FUN=function(x){strsplit(x,"")}))))

priors <- data.table(prior = strtrim(transitions[[1]], order), count = transitions[[2]])
priors <- aggregate(priors[[2]], by=list(priors[[1]]), FUN=sum)
priors <- data.table(word = priors[,1], count = priors[,2])

P.hat <- data.table(prior_state = strtrim(transitions[[1]], order),
    transition_state = substr(transitions[[1]], order+1, order+1),
    probability = apply(transitions, 1, 
      FUN=function(x){prior <- strtrim(x[[1]], order)
        return(as.numeric(x[[2]])/priors[word == prior,][[2]])})
)

initial <- table(apply(array(seqs), 1, FUN=function(y){strtrim(y,order)}))
initial <- data.table(word = names(initial), count = as.vector(initial))
initial <- initial[apply(array(initial[[1]]),1,nchar) == order]

n<-sum(initial[[2]])
alpha <- data.frame(alpha = (initial[[2]])/n, row.names = initial[[1]])

likelihood <- log(P.hat[[3]]) * transitions[[2]]
number_structures <- width(seqs)
sequences <- seqs[which(number_structures < order)]
number_structures <- width(sequences)
short_read_likelihood <- numeric(length=length(number_structures))
if(length(number_structures) > 0){
  for(x in 1:length(number_structures)){
    start <- sequences[x]
    n <- 0
    while(n == 0){
      state <- paste(start,"*",sep="")
      start <- strtrim(start,(nchar(start)-1))
      n <- sum(alpha[grepl(glob2rx(state), rownames(alpha)),])
    }
    short_read_likelihood[x] <- log(n*((1/length(structures))^(order-(nchar(state)-1))))
  }
}
likelihood <- sum(likelihood) + sum(log(alpha$alpha) * initial[[2]]) + sum(short_read_likelihood)

return(new("MC_Sequence_Analysis", 
          fasta_file=fasta_file,
          order=order,
          structures=as.factor(structures),
          alpha=alpha,
          transition_frequencies=P.hat,
          likelihood=likelihood
          )
      )
}
