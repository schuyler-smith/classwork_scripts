## BCB 568 HW #3

library(plyr)
library(gtools)
library(tidyverse)
library(reshape)
library(data.table)

setwd("C:/users/Schuyler/Dropbox/Classes/BCB_II/BCB_II_HW_3")
setwd("~/Dropbox/Classes/BCB_II/BCB_II_HW_3")

source("source.R")

### Reading data in 
d <- read.csv("count_genome.csv", header=TRUE)
d.l <- reshape(d, varying = paste("X", 2009:2016, sep=""), v.names = "counts", timevar = "yr", time = 2009:2016, direction = "long")
colnames(d.l) <- c("genome", "yr", "counts", "id")
d.ll <- untable(d.l, d.l$counts)[, -c(3,4)]

###1b
##### 2013:
d.l.2013 <- filter(d.l, yr==2013)  #choose only 2013
d.l.2013 <- arrange(d.l.2013, genome) %>% #make sure genomes are in order
  mutate(i = c(1:nrow(d.l.2013))) %>% #label data with the numbers "i"
  mutate(p.hat.i=counts/(sum(d.l.2013$counts)), #make a column of the MLE p.hat.i = xi/n
         var.p.hat.i=(p.hat.i)*(1-p.hat.i)/(sum(d.l.2013$counts)), 
         CI.90.upper = qnorm(0.95, mean=(p.hat.i), sd=(sqrt(var.p.hat.i))), #Make columns with 90% confidence interval bounds
         CI.90.lower = qnorm(0.05, mean=(p.hat.i), sd=(sqrt(var.p.hat.i))))

### plot MLE and confidence interval:
ggplot(d.l.2013, aes(i,p.hat.i)) +
  geom_pointrange(aes(ymin=CI.90.lower, ymax=CI.90.upper)) +
  ggtitle("2013")

### 2016:
d.l.2016 <- filter(d.l, yr==2016)  #choose only 2016
d.l.2016 <- arrange(d.l.2016, genome) %>% #make sure genomes are in order
  mutate(i = c(1:nrow(d.l.2016))) %>%
  mutate(p.hat.i=counts/(sum(d.l.2016$counts)), #make a column of the MLE p.hat.i = xi/n
         var.p.hat.i=(p.hat.i)*(1-p.hat.i)/(sum(d.l.2016$counts)),
         CI.90.upper = qnorm(0.95, mean=(p.hat.i), sd=(sqrt(var.p.hat.i))), #Make columns with 90% confidence interval bounds
         CI.90.lower = qnorm(0.05, mean=(p.hat.i), sd=(sqrt(var.p.hat.i))))

### plot MLE and confidence interval:
ggplot(d.l.2016, aes(i,p.hat.i)) +
  geom_pointrange(aes(ymin=CI.90.lower, ymax=CI.90.upper)) +
  ggtitle("2016")

### 1c

### Wald test statistic (T) for each year = sqrt((p_hat_i - p)^2/(var(p_hat_i)))
### Here, comparing genomes A and B (1 and 2, respectively)
(T.2013 <- sqrt(((d.l.2013$p.hat.i[1] - d.l.2013$p.hat.i[2])^2) / (d.l.2013$var.p.hat.i[1])))
(T.2016 <- sqrt(((d.l.2016$p.hat.i[1] - d.l.2016$p.hat.i[2])^2) / (d.l.2016$var.p.hat.i[1])))

# T ~ chi-square with 1 df, so find P(T<=t):
# (P.2013 <- pnorm(T.2013, 0, 1))
# (P.2016 <- pnorm(T.2016, 0, 1))
(P.2013 <- pchisq(T.2013, 1))
(P.2016 <- pchisq(T.2016, 1))

# Were A and B significantly different in either year? With two-tailed test, a=.05
(significant.2013 <- (P.2013 >= .975 | P.2013 < .025))
(significant.2016 <- (P.2016 >= .975 | P.2016 < .025))

# 2a #######
# Generate n bootstrap samples the same size as the original sample, with replacement

n <- 5
genome_letters <- LETTERS[1:length(levels(d.l$genome))]
results.2013 <- genome_letters
results.2016 <- genome_letters

for (i in 1:n){
  # Pick the 525 individuals for the resample.
  bootstrap_resample <- sample_n(d.ll,525, replace = TRUE)
  colnames(bootstrap_resample) <- c("genome","yr")
  
  # Make a table for 2013 data and calculate that year's p.i values
  bootstrap.2013 <- filter(bootstrap_resample, yr==2013)
  bootstrap.2013 <- create_bootstrap_data(bootstrap.2013)
  bootstrap.2013 <- mutate(bootstrap.2013, p=counts/sum(counts))
  p.2013 <- bootstrap.2013$p
  results.2013 <- cbind.data.frame(results.2013, p.2013)
  
  #Repeat above for 2016
  bootstrap.2016 <- filter(bootstrap_resample, yr==2016)
  bootstrap.2016 <- create_bootstrap_data(bootstrap.2016)
  bootstrap.2016 <- mutate(bootstrap.2016, p=counts/sum(counts))
  p.2016 <- bootstrap.2016$p
  results.2016 <- cbind.data.frame(results.2016, p.2016)
}

# Reorganize the results tables
names(results.2013) <- c("genome",1:n)
rownames(results.2013) <- results.2013$genome
results.2013 <- results.2013[,-1]
results.2013 <- t(results.2013)
results.2013 <- as.data.table(results.2013)

names(results.2016) <- c("genome",1:n)
rownames(results.2016) <- results.2016$genome
results.2016 <- results.2016[,-1]
results.2016 <- t(results.2016)
results.2016 <- as.data.table(results.2016)

# Find 90% confidence interval bounds for genotypes A and B, years 2013 and 2016

CI.90.2013 <- cbind.data.frame(genome = genome_letters, 
                lower.bound = apply(results.2013,2,FUN=function(x){quantile(x, .95)}),
                upper.bound =  apply(results.2013,2,FUN=function(x){quantile(x, .05)}),
                median =  apply(results.2013,2,FUN=function(x){quantile(x, .50)})
)

CI.90.2016 <- cbind.data.frame(genome = genome_letters, 
                               lower.bound = apply(results.2016,2,FUN=function(x){quantile(x, .95)}),
                               upper.bound =  apply(results.2016,2,FUN=function(x){quantile(x, .05)}),
                               median =  apply(results.2016,2,FUN=function(x){quantile(x, .50)})
)

# graph the confidence intervals to check for asymmetry
# 2013
ggplot(CI.90.2013, aes(i, median)) +
  geom_pointrange(aes(ymin=lower.bound, ymax=upper.bound)) +
  ggtitle("2013")

# 2016
ggplot(CI.90.2016, aes(i, median)) +
  geom_pointrange(aes(ymin=lower.bound, ymax=upper.bound)) +
  ggtitle("2016")

# 3b #############

# First sort genotypes by abundance.
sorted.2013 <- d.l.2013 %>% 
  select(genome, counts) %>% 
  arrange(desc(counts))

sorted.2016 <- d.l.2016 %>% 
  select(genome, counts) %>% 
  arrange(desc(counts))

grouping.info.2013 <- maximize.likelihood(sorted.2013, 3)
print(paste("Best grouping for 2013 had likelihood", grouping.info.2013[1]))
print_groups(sorted.2013,grouping.info.2013)

grouping.info.2016 <- maximize.likelihood(sorted.2016, 3)
print(paste("Best grouping for 2013 had likelihood", grouping.info.2016[1]))
print_groups(sorted.2016,grouping.info.2016)

# The above stuff was done with genotypes with 0 abundance.  From an experimental standpoint,
# we should include the genotypes with 0 observations because we do want to make inferences about
# their abundance and put them into groups.  However, since our log-likelihood equation
# takes the log of each lambda, the values of the log-likelihood will tank when lambda = 0 and
# that will mess up the way we choose our groups.  Below is the above calculations done by removing
# genotypes with 0 observations from the 2013 and 2016 data.

sorted.2013.trimmed <- sorted.2013 %>% filter(counts != 0)
sorted.2016.trimmed <- sorted.2016 %>% filter(counts != 0)

grouping.info.2013.trimmed <- maximize.likelihood(sorted.2013.trimmed, 3)
print(paste("Best grouping for 2013 had likelihood", grouping.info.2013.trimmed[1]))
print_groups(sorted.2013.trimmed,grouping.info.2013)

grouping.info.2016.trimmed <- maximize.likelihood(sorted.2016.trimmed, 3)
print(paste("Best grouping for 2016 had likelihood", grouping.info.2016.trimmed[1]))
print_groups(sorted.2016.trimmed,grouping.info.2016)
#These groups seem more reasonable.
