options(warn=-1)
# local.libpath <- "C:/Program Files/R/R-3.4.3/library"
local.libpath <- "~/R/x86_64-pc-linux-gnu-library/3.4"
.cran_packages <- c("devtools", "curl","dplyr","reshape2", "tidyverse", "magrittr", "readr", "stringr", "data.table", "svMisc")
.bioc_packages <- c("Biostrings")
.local_packages <-c("zwc")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}; sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
.inst <- .local_packages %in% installed.packages()
if(any(!.inst)) {
  withr::with_libpaths(new = local.libpath, devtools::install_github('arendsee/zwc'))
}; sapply(c(.local_packages), require, character.only = TRUE)

setwd("~/Dropbox/Classes/BCB_II/BCB_II_HW_4")
# setwd("C:/Users/schuyler/Dropbox/Classes/BCB_II/BCB_II_HW_4")
source("Functions/MC_seq_analysis.R")
source("Functions/simulate_MC_reads.R")

# Set the dataset for analysis, either test or real data
#test pdb test-v2 pdb-v2

fasta <- "pdb"
order <- 3

fasta_file <- paste("Fastas/",fasta,".fa",sep="")

order_analysis <- MC_seq_analysis(fasta_file = fasta_file, order = order)
saveRDS(order_analysis, file = paste("MC_Analysis/order_",order,"_",fasta,".RDS", sep=""))


# simulated_reads <- simulate_MC_reads(order_analysis)
# writeXStringSet(BStringSet(simulated_reads), "Fastas/simulated_reads_order_2.fa")

