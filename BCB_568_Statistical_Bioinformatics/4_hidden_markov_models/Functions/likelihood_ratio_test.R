
likelihood_ratio_test<- function(likelihood_test, confidence = 0.95){
  likelihood_test[, df := c(NA, apply(array(1:(length(likelihood_test[[2]])-1)), 1, FUN=function(i){7*(8^(i+1))}))]
  likelihood_test[, t := c(NA, round(2*diff(likelihood_test[[2]])))]
  likelihood_test[, chisqr := c(NA, apply(array(1:(length(likelihood_test[[2]])-1)), 1, FUN=function(i){round(qchisq(p=confidence, df = 7*(8^(i+1))))}))]
  return(likelihood_test[, reject_H0 := array(apply(likelihood_test,1,FUN=function(i){i[[4]] > i [[5]]}))])
}

