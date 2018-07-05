library(tidyverse)
library(plyr)
library(dplyr)

setwd("~/Dropbox/Classes/BCB_II/BCB_II_HW_2")
source("functions.R")

seq_data <- read.table("ngs_site_data.Rtxt", header = T)
seq_data <- mutate(seq_data, "Pg0" = case_when(a == 1 ~ prob_match(q, 0), a == 0 ~ prob_mismatch(q, 0), TRUE ~ 0),
                   "Pg1" = case_when(a == 1 ~ prob_match(q, 1), a == 0 ~ prob_mismatch(q, 1), TRUE ~ 0),
                   "Pg2" = case_when(a == 1 ~ prob_match(q, 2), a == 0 ~ prob_mismatch(q, 2), TRUE ~ 0))

# Multiply all the conditional match|mismatch probabilities like 3c. 
cond_probs = ddply(seq_data, ~i, summarize, p0 = prod(Pg0), p1 = prod(Pg1), p2 = prod(Pg2))
# 8a
plot(Vectorize(log.likelihood), xlab = "psi", ylab = "log-likelihood(psi)", ylim = c(-120, -40))
# 8b
psihat <- optim(0.5, lower = 0, upper = 1, log.likelihood, method = 'Brent', control = list(fnscale = - 1))$par

#9
bootstrap_psi(1, cond_probs)
#10
bootstrapped_psi_estimates <- bootstrap_psi(10000, cond_probs)
hist(bootstrapped_psi_estimates, xlab = "psi", main = "Frequency of psi estimates from bootstrapping")
# Getting a confidence interval, we can approximate a 90% confidence interval by sorting our list of psi estimates
# and picking the 501st least and 950th greatest elements.  
sorted_bootstrap_estimates <- sort(bootstrapped_psi_estimates)
CI_90 <- c(sorted_bootstrap_estimates[501], sorted_bootstrap_estimates[9500])
CI_90

# 11
confidence_level <- 0.90 
# Find corresponding eta_p, the pth quantile of the chi square distribution
eta <- qchisq(confidence_level, df=1) 
lower_bound <- uniroot(lambda, c(.5,.8), tol= 0.0001)$root
upper_bound <- uniroot(lambda, c(.8, .9), tol= 0.0001)$root
lower_bound
upper_bound










