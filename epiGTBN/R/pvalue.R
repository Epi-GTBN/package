pvalue <- function(mylist1,mylist2,mylist3) {
	mylist1 <- data.frame(mylist1)
	mylist2 <- data.frame(mylist2)
	mylist3 <- data.frame(mylist3)
	df1 <- nrow(unique(mylist1))
	df2 <- nrow(unique(mylist2))
	df3 <- nrow(unique(mylist3))
	df <- (df1-1) * (df2-1) * df3
	cmi <- entropy(cbind(mylist1,mylist3)) + entropy(cbind(mylist2,mylist3)) - entropy(mylist3) - entropy(cbind(mylist1,mylist2,mylist3))
	return(pchisq(2 * cmi * nrow(mydata),df,FALSE,FALSE))
}