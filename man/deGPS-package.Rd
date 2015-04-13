\name{deGPS-package}
\alias{deGPS-package}
\alias{deGPS}
\docType{package}
\title{
Normalization and Two-group Differential Expression Test for RNA-seq Data
}
\description{
This package is proposed to analyze RNA-seq data in two steps: normalization and differential expression test. 

New normalization methods based on generalized Poisson distribution are contained in the main analysis functions, in which \code{GP-Theta} is suggested. Other popular normalization methods, such as TMM, LOWESS, Quantile, is also availbale in the package. More than one method can be specified in one run, in which case the resulting p values are a p value matrix, with each column representing one method.

Differential expression test is designed for RNA-seq data, especially those of small sample size, based on permutation strategy. Note that deGPS can only handle DE tests between two groups, but the novel GP-based normalization methods can be applied on any RNA-seq data. More Details about the GP-based normalization method and the advantages and limitations of our DE test can be found in our article.

The permutation step may become time-consuming in large sample size context or in mRNA read count data analysis. Parallel computation is introduced in the main functions to deal with the computational burden. Note that in situations where parallel computation is not necessary, to force parallelling may result in more run time.

The package also contains function to generate GP distrbuted data to be an example of RNA-seq data. It has to be pointed out that the simulated data in this package is FAR AWAY from real RNA-seq data, it is just simple GP distributed samples. For more appropriately simulated RNA-seq, \code{compcodeR} package is suggested, or, alternatively, you can generate H0 and H1 data from real data. A simple example is given in the following example session. More R codes referring to the real data based simulation or compcodeR based simulation can be found in the supplementary materials of our article.

Please do not hesitate to contact the author if you have any questions or find any flaws in the package.
}
\details{
\tabular{ll}{
Package: \tab deGPS\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-02-24\cr
License: \tab GPL-2\cr
}
To do DE test for miRNA, call \link{GPSmle}. To do DE test for mRNA, call \link{deGPS_mRNA}. If only the normalization step is required, call \link{GPSmle} and specify \code{type} as "normalization". 

Version 2.0 includes \code{sam.stat} in package \code{st} in \link{deGPS_mRNA} to get empirical t distribution instead of traditional t stats.

Some denotations used in this manual are as follow:
\itemize{
  \item \code{nSample} 

sample size of the data. If there are more than one replicate for each sample, nSample represents the sample size multiplied by the number of replicates.
  \item \code{nRNA} 

the number of genes in the data.
  \item \code{n-miRNA} 

the number of miRNAs in the data.
  
  \item \code{groupSize} 

the number of samples in either group.
  \item \code{permutationTimes} 

the times of permutation used to obtain the empirical distribution of T-stat.

  \item \code{round(.)} 

the round function in R. The value of \code{round(x)} is the closest integer of x.
  \item \code{ceiling(.)} 

the ceiling function in R. The value of \code{ceiling(x)} is the integer part of x added by 1.
  \item \code{ncol(.)} 

the ncol function in R. The value of \code{ncol(data)} is the number of columns of data.
  \item \code{mod(.)} 

the mod function in R. The value of \code{mod(x)} is the remainder of x.
}
}
\author{
Chen Chu

Maintainer: Chen Chu <chuchenblueblues@gmail.com>
}
\references{
deGPS: a Powerful and Flexible Framework for Detecting Differential Expression in RNA-Sequencing Studies
}
\keyword{ GP }
\keyword{ RNA-seq }
\keyword{ DifferentialExpression}
\examples{
\dontrun{
### Analyze real RNA-seq data
### Compare "Early Embryo" to "Late Embryo" in fly data
data(flyData)

i <- 1

group <- rep(1:2, each = 6)
simuData <- as.matrix(flyData$data[ , -1])
simuData <- simuData[ , flyData$compIdx[[i]]]


###remove genes of all-zero read counts
simuData <- simuData[apply(simuData, 1, function(x) !all(x == 0)), ]    

### apply deGPS on fly data with the empirical T stats downloaded at 
### https://www.dropbox.com/s/if5ido5vd8rzff5/empirical.T.stats.fly.RData
### you can set empirical.T.stats as NULL to get your own empirical values
### note that for non-parallelized computation, it may take hours.

load("empirical.T.stats.fly.RData")
empirical.T.stats.fly <- empTAll[i]
names(empirical.T.stats.fly) <- "GP-Theta"

### Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
flyRes <- deGPS_mRNA(data = simuData, group = group, 
ncore = 6, nSubcore = 6, method = "GP-Theta", 
empirical.T.stats = empirical.T.stats.fly)

pvalue <- as.vector(flyRes$pvalue)
adj.p <- p.adjust(pvalue, method = "BH")
topTags(pvalue)

### Generate Random samples from GP(theta, lambda)
examData <- newExampleData(nRNA = 100, groupSize = 6, lambda = 0.9, 
theta = 3, ptol = 1e-15)
str(examData)

### Differential Expression Tests
examRes <- deGPS_mRNA(data = examData$data, group = examData$group, 
method = "GP-Theta", nSubcore = 2, ncore = 2, geneid = paste("G", 1:100, sep = ""))
str(examRes)

### Generate simulated RNA-seq data from compcodeR package
require(compcodeR)

samples.per.cond <- 5
random.outlier.high.prob <- 0.1
n.vars <- 10000

examData <- generateSyntheticData(dataset = "simuData",
n.vars = n.vars, samples.per.cond = samples.per.cond, n.diffexp = floor(n.vars * 0.1),
repl.id = 1, seqdepth = 1e+07, fraction.upregulated = 0.5,
between.group.diffdisp = FALSE, filter.threshold.total = 1,
filter.threshold.mediancpm = 0, fraction.non.overdispersed = 0,
random.outlier.high.prob = random.outlier.high.prob,
output.file = "simuData_repl1.rds")

group <- examData@sample.annotations$condition

### Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
examRes <- deGPS_mRNA(data = examData@count.matrix, group = group, 
method = "GP-Theta", nSubcore = 6, ncore = 6, geneid = paste("G", 1:nrow(examData@count.matrix), sep = ""))
str(examRes)

}
}