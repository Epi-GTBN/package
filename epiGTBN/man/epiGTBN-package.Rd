\name{epiGTBN-package}
\alias{epiGTBN-package}
\alias{epiGTBN}
\docType{package}
\title{An Approach of Epistasis Mining Based on Genetic Tabu Algorithm and Bayesian Network}
\description{

  Learn the structure of a Bayesian network using genetic tabu algorithm for 2-locus epistasis and
  3-locus epistasis.

}
\details{

  \tabular{ll}{

    Package: \tab epiGTBN\cr
    Type: \tab Package\cr
    Version: \tab 2.1.0-3\cr
    Date: \tab 2019-12-11\cr
    License: \tab GPLv2 or later

  }

  This package implements an approach of epistasis mining based on Genetic Tabu Algorithm and Bayesian Network.
  It uses Genetic Tabu Algorithm into the heuristic search strategy of Bayesian network.

  The individual structure can be evolved through the genetic manipulations of selection, crossover and mutation.
  It can help to find the optimal network structure, and then mine the epistasis loci in further.

  In order to enhance the diversity of the population and obtain a more effective global optimal solution, we implement
  the tabu search strategy into the operations of crossover and mutation in genetic algorithm.

  For additional information and supplementary material, visit http://122.205.95.139/Epi-GTBN/.
}

\author{

  Yang (Simon) Guo, ZhiMan (Melody) Zhong\cr
  College of Informatics\cr
  Huazhong Agricultural University (HZAU)

  Maintainer: Yang (Simon) Guo \email{simon-guo@outlook.com}

}

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
\keyword{ package }
