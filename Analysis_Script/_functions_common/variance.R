variance <- function(x) {
  
  val <- var(x)*(length(x)-1)/length(x)
  return(val)
  
}