# setwd("")
source("Functions/likelihood_ratio_test.R")
orders <- 8
fasta <- "simulated_reads_order_2"


likelihood_table <- data.table(order = 1:orders,
            LL = apply(array(1:orders), 1, FUN=function(x){
                sum(readRDS(paste("MC_Analysis/order_",x,"_",fasta,".RDS", sep=""))@likelihood)}))
likelihood_ratio_test(likelihood_table, 0.95)
likelihood_table


likelihood_table <- data.table(read.table("MC_Analysis/pdb.likelihoods.trimmed", header=TRUE ,sep=","))
likelihood_table <- likelihood_table[-c(9,10),]
likelihood_ratio_test(likelihood_table, 0.95)
colnames(likelihood_table) <- c("order", "LL","df","t","chisqr","reject_H0")
likelihood_table

