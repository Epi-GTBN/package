# binary mutual information calculation frontend.
# for 2 nodes

# dev2 branch, has to convert R numeric matrix to C 0/1 array in every arc MI calc, time-consuming.
# extract convertion part to a new function `num2binary` in dev3 branch
# mi.2 <- function(X, varNum1, varNum2, varNum3) {
#     # print(X)
#     nsamples <- NCOL(X) 
#     nvars <- NROW(X) 
#     # cat("nsamples: ", nsamples, " nvars: ", nvars, "\n")
#     res <- NULL
#     res <- .Call("mi2R", X, varNum1, varNum2, varNum3, nsamples, nvars)
#     res
# }

# dev3@branch: frontend for R numeric matrix to C 0/1 array transformation
mi.2 <- function(X, varNum2, varNum3) {
    # print(X)
    nsamples <- NCOL(X)
    # cat("nsamples: ", nsamples, "\n")
    res <- NULL
    res <- .Call("mi2R", varNum2, varNum3, nsamples)
    res
}