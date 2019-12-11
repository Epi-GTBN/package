
entropy <- function( X, method = "emp")
{
      X<-data.frame(X)
	  X<-data.matrix(X)
	  n <- NCOL(X) #nsamples
      N <- NROW(X) #nvars
	  Z<-na.omit(X)
	  if( !(all(Z==round(Z))))
	      stop("This method requires discrete values")                      
	#if(n>32000)
		#stop("too many variables")
    res <- NULL 

    if( method == "emp")
		choi<-0
	else if( method == "mm" )
		choi<-1
	else if( method == "sg" )
		choi<-2
	else if(method == "shrink")
		choi<-3
	else stop("unknown method")
	res <- .Call( "entropyR",X,N,n, choi)
    # SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice
    #             int nsamples, int nvars
	res
}
