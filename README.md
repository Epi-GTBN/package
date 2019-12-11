<!-- <img align="center" src='./HZAU-120th-slogan.jpg'/> -->
<img align="center" src='./Epi-GTBN_banner.png'/>

<h1 align="center">
  R package for Epi-GTBN
</h1>

<h4 align="center">An Approach of Epistasis Mining Based on Genetic Tabu Algorithm and Bayesian Network</h4>

<p align="center">
  Current version of this package: <strong>2.1.0-2</strong><br/>
  <a href="https://cran.r-project.org"><img src="https://img.shields.io/badge/R-3.4.1-green.svg" alt="R-3.4.1"/></a>
  <a href="http://gcc.gnu.org"><img src="https://img.shields.io/badge/gcc-7.1.0-green.svg" alt="gcc-7.1.0"/></a>
</p>

This package implements an approach of epistasis mining based on Genetic Tabu Algorithm and Bayesian Network.
It uses Genetic Tabu Algorithm into the heuristic search strategy of Bayesian network.

The individual structure can be evolved through the genetic manipulations of selection, crossover and mutation. It can help to find the optimal network structure, and then mine the epistasis loci in further.

In order to enhance the diversity of the population and obtain a more effective global optimal solution, we implement the tabu search strategy into the operations of crossover and mutation in genetic algorithm.                    

## INSTALLATION
                    
- For Mac and Linux user: enter directory that stores the `.tar.gz` package, then run command `R CMD install epiGTBN` to install
- For PC user: additional software 'Rtools' needs to be installed, then in R console run command `install.packages("path\\to\\epiGTBN_XXX.tar.gz", repos=NULL, type="source")` to install
- Install from github: TBA

## USAGE EXAMPLE

You can also use `help("epiGTBN-package")`, `help("epiGTBN")`, `help("gtbn2")`, `help("gtbn3")` in R to access usage example.

```R
library(epiGTBN)
# ------------------------------
# For 2-locus epistasis mining, use function gtbn2(). e.g.
# load testing data
data(Epistasis_Example_2_locus)
# load parameters for GTBN
# note that `parameters.R` should be placed in working diretory, use getwd() or setwd() to get or set your working diretory
source("parameters.R")
# max.iter indicates max iteration before stopping epiGTBN
res1 <- gtbn2(Epistasis_Example_2_locus, max.iter = 60, debug = FALSE)
# epistasis mining result is presented as arcs in graph, it can be export to data frame using code below
prediction <- data.frame(res1$arcs)
# ------------------------------

# For 3-locus epistasis mining, use function gtbn3() instead. e.g.
# ------------------------------
# load testing data
data("Epistasis_Example_3_locus")
# load parameters for GTBN
# note that `parameters.R` should be placed in working diretory, use getwd() or setwd() to get or set your working diretory
source("parameters.R")
# max.iter indicates max iteration before stopping epiGTBN
res2 <- gtbn3(Epistasis_Example_3_locus, max.iter = 60, debug = FALSE)
# epistasis mining result is presented as arcs in graph, it can be export to data frame using code below
prediction <- data.frame(res2$arcs)
# ------------------------------

# If you wish to see debug output, set debug = TRUE and use sink(). e.g. 
# ------------------------------
sink(file = "GTBN-2nodes-log.txt")
res1 <- gtbn2(Epistasis_Example_2_locus, max.iter = 60, debug = TRUE)
sink()
# and/or
sink(file = "GTBN-3nodes-log.txt")
res2 <- gtbn3(Epistasis_Example_3_locus, max.iter = 60, debug = TRUE)
sink()
# ------------------------------

# Blacklist and whitelist use
# ------------------------------
# e.g. the arc N1 - N8 should not be there?
blacklist <- data.frame(from = c("N1", "N8"), to = c("N8", "N1"))
blacklist
res3 <- gtbn2(..., blacklist = blacklist, ...)

# e.g. force N1 - N9 direction (N1 -> N9).
whitelist <- data.frame(from = c("N1"), to = c("N9"))
whitelist
res4 <- gtbn2(..., whitelist = whitelist, ...)

# e.g. use both blacklist and whitelist.
res5 <- gtbn2(..., whitelist = whitelist, blacklist = blacklist, ...)
```

## INPUT FORMAT                    

The last column of input data should be phenotype Class, where 1 represents case, 0 represents control. 
In the rest columns, (0, 1, 2) is used to express the genotype data. Specifically, 0 denotes homozygote common genotype, 1 denotes heterozygous genotype and 2 denotes homozygote rare genotype.

## LICENSING

Epi-GTBN is distributed under GPL version 2 or later, see the source code.    

## Publication

> Yang Guo#, Zhiman Zhong#, Chen Yang, Jiangfeng Hu, Yaling Jiang, Zizhen Liang, Hui Gao, Jianxiao Liu*. Epi-GTBN: An Approach of Epistasis Mining Based on Genetic Tabu Algorithm and Bayesian Network. BMC Bioinformatics 20(1): 444:1-444:18 (2019)

[Paper (.pdf)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3022-z)

## Organization

<img src="./HZAU-120th.png" alt="HZAU-120th" width = 340px/><br/>

<img src="./COI.png" alt="College of Infomatics, HZAU" />

<a href="https://github.com/Epi-GTBN"><img src="https://sgyzetrov.github.io/images/epiGTBN-horizontal.png" alt="Epi-GTBN logo" title="An Approach of Epistasis Mining Based on Genetic Tabu Algorithm and Bayesian Network" height = 100px/></a>

