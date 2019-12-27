mi_partial <- function(CachedMatrix3Col,CachedMatrix4Col,mylist2,mylist3,mylist4){ 
  # Two methods to improve the speed of calculating joint mutual information in a triple loop:
  # 1. Calculate entropy (mydata [, 101]) outside the loop, reducing the time required to calculate entropy (mydata [, 101]) every time the conditions are met
  # 2. No longer use time-consuming cbind, but rather allocate memory for objects in advance
  CachedMatrix3Col[,1] <- mylist2
  CachedMatrix3Col[,2] <- mylist3
  CachedMatrix3Col[,3] <- mylist4
  CachedMatrix4Col[,c(2,3,4)] <- CachedMatrix3Col
  return(entropy(CachedMatrix3Col) - entropy(CachedMatrix4Col))
}