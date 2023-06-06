# calculate effect size

EffectSize <- function(m1, sd1, m2, sd2) {
  
  res <- (m1 - m2)/sqrt((sd1^2 + sd2^2)/2)
  return(res)
  
}