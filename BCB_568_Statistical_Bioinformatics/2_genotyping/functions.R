phred_to_error <- function(Q){return(10^(-Q/10))}

genotype_prob <- function(psi, g){ 
  return(dbinom(g, 2, psi))
}

prob_match <- function(Q, g){
  e <- phred_to_error(Q)
  match <- e*(2-g)/3 
  mismatch <- g*(1-e)
  return(round(((match + mismatch)/2),2))
}
prob_mismatch <- function(Q, g){
  e <- phred_to_error(Q)
  mismatch <- e*g/3 
  match <- (2-g)*(1-e)
  return(round(((mismatch + match)/2),2))
}

log.likelihood <- function(psi, conditional_probabilities = cond_probs){
  return(sum(apply(conditional_probabilities[,2:4], 1, FUN=function(x){log(sum(x*genotype_prob(psi, 0:2)))})))
}

bootstrap_psi <- function(n = 1, samp = cond_probs){
  return(sapply(vector("numeric", length = n), FUN=function(x){optim(0.5, lower = 0, upper = 1, log.likelihood, method = 'Brent', conditional_probabilities = samp[sample(10, replace = TRUE), ], control=list(fnscale = -0.01))$par
  }))
}

lambda <- function(psi0){
  -2*(log.likelihood(psi0) - log.likelihood(psihat)) - eta 
}

