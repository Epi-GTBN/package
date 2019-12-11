mi_partial <- function(CachedMatrix3Col,CachedMatrix4Col,mylist2,mylist3,mylist4){ 
  # 两种方式提高了3重循环中计算联合互信息的速度：
  # 1. 在循环外计算entropy(mydata[,101])，减少每次符合条件时都需要算一遍entropy(mydata[,101])的时间
  # 2. 不再使用十分耗时的cbind，而是预先给对象分配内存
  CachedMatrix3Col[,1] <- mylist2
  CachedMatrix3Col[,2] <- mylist3
  CachedMatrix3Col[,3] <- mylist4
  CachedMatrix4Col[,c(2,3,4)] <- CachedMatrix3Col
  return(entropy(CachedMatrix3Col) - entropy(CachedMatrix4Col))
}