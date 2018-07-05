log.likelihood <- function(lambda, x)
{
  # if lambda=0, then log(lambda) = -Inf, which will cause calculation problems.
  # Instead, set log.likelihood to a large negative number if lambda==0.
  # See discussion on assignment.
  if (lambda == 0){
    return(-100000)
  }
  return(x * log(lambda) - lambda - log(factorial(x)))  
}

maximize.likelihood <- function(data_set, K){
  partitions <- combinations(nrow(data_set) - 1, K - 1)
  best.likelihood = -1000000
  
  for (i in 1:nrow(partitions))
  {
    breaks <- c(0,partitions[i,],nrow(data_set))
    groups <- lapply(1:K, FUN=function(x){
      data_set[(breaks[x]+1):breaks[x+1], 2]
    })
    group_members <- c(apply(array(1:K), 1, FUN=function(x){
      length(data_set[(breaks[x]+1):breaks[x+1], 1])
    }))
    lambdas = sapply(groups, sum)
    current.likelihood <- sum(sapply(1:K, FUN = function(x){
      sum(log.likelihood(lambdas[x], groups[[x]]))
    }))
    
    if (current.likelihood > best.likelihood) {
      best.likelihood = current.likelihood
      best.partition = breaks[-c(1,length(breaks))]
      best.lambdas = lambdas
      best.group.members = group_members
    }
  }
  return(list(Max_Likelihood = best.likelihood, Partitions = best.partition, Lambdas = best.lambdas, Group_Members = best.group.members))
}

print_groups <- function(data_set, group_info){
  print("Groups: ")
  breaks = c(0, group_info$Partitions, nrow(data_set))
  groups <- apply(array(1:length(group_info$Lambdas)),1,FUN=function(x){
    print(data_set[(breaks[x]+1):breaks[x+1], 1], max.levels = 0)
  })
}

create_bootstrap_data <- function(data_set){cbind.data.frame(genome = LETTERS[1:19], counts = apply(array(LETTERS[1:length(levels(data_set$genome))]), 1, FUN=function(x){
  sum(data_set$genome==x)
}))}
