options(warn=-1)
.cran_packages <- c("kernlab","smotefamily", "pROC", "stringr", "data.table", "FNN")
.bioc_packages <- c("")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {install.packages(.cran_packages[!.inst])}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {source("http://bioconductor.org/biocLite.R");biocLite(.bioc_packages[!.inst], ask = F)}; sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
