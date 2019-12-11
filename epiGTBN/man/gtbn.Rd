\name{score-based algorithms}
\alias{score-based algorithms}
\alias{gtbn2}
\alias{gtbn3}
\title{An Approach of Epistasis Mining Based on Genetic Tabu Algorithm and Bayesian Network}
\description{

  Learn the structure of a Bayesian network using genetic tabu algorithm for 2-locus epistasis (gtbn2) or
  3-locus epistasis (gtbn3).
  
  For additional information and supplementary material, visit http://122.205.95.139/Epi-GTBN/.
}
\usage{
gtbn2(x, start = NULL, whitelist = NULL, blacklist = NULL, max.iter = Inf,
  debug = FALSE)
gtbn3(x, start = NULL, whitelist = NULL, blacklist = NULL, max.iter = Inf,
  debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{start}{an object of class \code{bn}, the preseeded directed acyclic
    graph used to initialize the algorithm. If none is specified, an empty one
    (i.e. without any arc) is used.}
  \item{whitelist}{a data frame with two columns (optionally labeled "from" and
    "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labeled "from" and
    "to"), containing a set of arcs not to be included in the graph.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{max.iter}{an integer, the maximum number of iterations.}
}
\value{

  An object of class \code{bn}.

}
\author{ Yang (Simon) Guo }
\examples{
library(epiGTBN)
# load testing data
data(Epistasis_Example_2_locus)
# load parameters for GTBN
# note that <parameters.R> should be placed in working diretory, use getwd() or setwd() to get or set your working diretory
source("parameters.R")
# max.iter indicates max iteration before stopping GTBN
res1 <- gtbn2(Epistasis_Example_2_locus, max.iter = 60, debug = FALSE)
# epistasis mining result is presented as arcs in graph, it can be export to data frame using code below
prediction <- data.frame(res1$arcs)
#========================================================
# For 3-locus epistasis mining, use function gtbn3 instead
# load testing data
data("Epistasis_Example_3_locus")
# load parameters for GTBN
# note that <parameters.R> should be placed in working diretory, use getwd() or setwd() to get or set your working diretory
source("parameters.R")
# max.iter indicates max iteration before stopping GTBN
res2 <- gtbn3(Epistasis_Example_3_locus, max.iter = 60, debug = FALSE)
# epistasis mining result is presented as arcs in graph, it can be export to data frame using code below
prediction <- data.frame(res2$arcs)

# If you wish to see debug output, set debug = TRUE and use sink()
\dontrun{
sink(file = "GTBN-2nodes-log.txt")
res1 <- gtbn2(Epistasis_Example_2_locus, max.iter = 60, debug = TRUE)
sink()
#========================================================
sink(file = "GTBN-3nodes-log.txt")
res2 <- gtbn3(Epistasis_Example_3_locus, max.iter = 60, debug = TRUE)
sink()
}

## Blacklist and whitelist use
# the arc N1 - N8 should not be there?
blacklist <- data.frame(from = c("N1", "N8"), to = c("N8", "N1"))
blacklist
res3 <- gtbn2(..., blacklist = blacklist, ...)
# force N1 - N9 direction (N1 -> N9).
whitelist <- data.frame(from = c("N1"), to = c("N9"))
whitelist
res4 <- gtbn2(..., whitelist = whitelist, ...)
# use both blacklist and whitelist.
res5 <- gtbn2(..., whitelist = whitelist, blacklist = blacklist, ...)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
