# calculate effect size
# Cohen's d when sample size of two groups is equal

EffectSize <- function(m1, sd1, m2, sd2) {
  
  res <- (m1 - m2)/sqrt((sd1^2 + sd2^2)/2)
  return(res)
  
}

EffectSize_n <- function(m1, sd1, m2, sd2, n1, n2) {
  
  res <- (m1 - m2)/sqrt((n1 * sd1^2 + n2 * sd2^2)/(n1 + n2))
  return(res)
  
}