
EM_trim_calc <- function(
  M = fasta_file,
  structures = c("A","C","G","T"),
  adapter = "GCCTTGCCACACGCTCAG",
  barcode = 'NNNNNNNNN',
  primer = 'GTTGTAAYTTCTAGRTCCCCTCCTG',
  Z = 3,
  q_Sb = rep(0.25, length(structures)),
  q_Bb = rep(0.25, length(structures)),
  gamma = rep(0.33, (length(structures)*(length(structures)))),
  delta = rep(0.99, sum(nchar(c(adapter,barcode,primer)))),
  xi_1 = rep(0.2, (Z+2)),
  xi_2 = rep(0.25, length(structures)),
  error = 0.01
){
  K = Z+2
  reference <- unlist(strsplit(paste(adapter, barcode, primer, sep=""),""))
  biallelic_regions = which((unlist(strsplit(reference,"")) %in% c("Y","R")))
  N_regions = which((unlist(strsplit(reference,""))) == "N")
  IP <- which(reference %in% structures)
  l_val = length(biallelic_regions)^2
  trim_lengths <- data.table(read = names(readBStringSet(M)))
  trim_lengths[, trim := rep(0, length(trim_lengths[[1]]))]
  trim_lengths[, k := rep(4, length(trim_lengths[[1]]))]
  trim_lengths[, E := rep(1, length(trim_lengths[[1]]))]
  EM_reads <- paste(readBStringSet(M))[width(paste(readBStringSet(M))) >= length(reference)+(K-1)]

  test <- function(x,y){
    return(((x-y)/(x+(error/1000))) <= error)}
  
  E <- matrix(0,length(EM_reads),((K-1)*l_val)+1)
  current <- 10
  previous <- 1
  loop <- 0
  while(any(mapply(test, current, previous) == FALSE)){
      previous <- E
      E <- calc_E(M = EM_reads, reference, K, l_val, q_Sb, q_Bb, gamma, delta, xi_1, xi_2, structures, IP)
      q_Sb <- calc_qSb(EM_reads, K, l_val, structures, E, reference)
      q_Bb <- calc_qBb(EM_reads, K, l_val, structures, E, N_regions)
      gamma <- calc_gamma(EM_reads, structures, reference, biallelic_regions, gamma, E, K, l_val, IP)
      delta <- calc_delta(EM_reads, reference, K, l_val, IP, E, biallelic_regions)
      xi_1 <- calc_xi1(EM_reads, E, K, l_val)
      xi_2 <- calc_xi2(E, K, l_val)
      current <- E
      loop <- loop+1
      print(paste(loop, " iterations completed", sep=""))
  }
  trim <- apply(E,1,FUN=function(x){which(x == max(x))})
  trim <- sapply(trim,FUN=function(x){floor((x-0.5)/4)})
  max_E <- apply(E,1,max)
  trim_lengths[which(width(paste(readBStringSet(M))) >= length(reference)+(K-1))][[2]] <- (trim + length(reference))
  trim_lengths[which(width(paste(readBStringSet(M))) >= length(reference)+(K-1))][[3]] <- trim
  trim_lengths[which(width(paste(readBStringSet(M))) >= length(reference)+(K-1))][[4]] <- max_E
  trim_lengths[k == 4][[2]] <- 0
  return(list(MLE_trim = trim_lengths, q = data.table( q_Sb = q_Sb, q_Bb = q_Bb), gamma = gamma, delta = delta))
}
