mi<-function(mylist1,mylist2,mylist3){
  return(entropy(mylist1)+entropy(cbind(mylist2,mylist3))-entropy(cbind(mylist1,mylist2,mylist3)))
}