\name{GPSmle.default}
\alias{GPSmle.default}
\title{
Generalized Poisson Statistical Maximum Likelihood Estimation (default)
}
\description{
the default method for the function \code{GPSmle}.
}
\usage{
\method{GPSmle}{default}(data, group = rep(1:2, each = 5), 
type = c("pvalue", "normalization", "ecdf"), 
method = c("GP-Theta", "Lowess", "GP-Quantile", "Quantile", "TMM"), 
maxIter = 500, paired = FALSE, ncpu = 1, geneid = NULL, empirical.T.stats = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 
a matrix containing microRNA read count data. The column represents samples while the row represents miRNAs. Biologcial or technical replicates must be made as new columns in the data.
}
  \item{group}{
The group indicator. The length of group must equal to ncol(data).
}
  \item{type}{
\code{type} can be "normalization", "ecdf" or "pvalue", to which step \code{GPSmle} will stop. "normalization" means that only normalized data is returned and no DE test is conducted; "ecdf" means the empirical T-stats are generated after nomalization, and the output contains both the normalized data sets and empirical values; "pvalue" means p-values are calculated after the empirical T-stats are obtained, and the output contains the normalized data sets, empirical T-stats and the p-values of DE test. 
}
  \item{method}{
The methods of normalization. More than one method can be specified. The value can be "Lowess", "Quantile", "TMM", "GP-Quantile", "GP-Theta" or "GP-MLE2L". See the reference for more details.
}
  \item{maxIter}{
The default value of maxIter is 500. When sample size is large, instead of transversing every possible permutations, randomly sampling is applied for \code{maxIter} to obtain the empirical distributions. Larger \code{maxIter} costs longer run time. Note that \code{maxIter} is forced to be not larger than \code{permutationTimes}.
}
  \item{paired}{
The current version of deGPS only contain unpaired test.
}
  \item{ncpu}{
The number of cores for the parallel computing. When sample size is large, the permutation step may be time-consuming. Specify \code{ncpu} > 1, parallel computation is applied in the function. \code{ncpu} cores are used to calculate \code{maxIter} times of permutations, each of which take responsibilies of part of the permutation task. The calculation of p values are also splitted into subsets for parallel computation.
}
  \item{geneid}{
Gene id of the specified data. Biological or technical replicates must be new columns in the data, i.e., duplicates in geneid are not allowed.
}
  \item{empirical.T.stats}{
A list of empirical T statistics with names as normalization methods, the format of which must be the same as the \code{empirical.T.stats} in the returned list of \code{deGPS_mRNA} or \code{GPSmle}, i.e., a list of empirical T-stats with the method names as the list name. If null, the empirical T statistics will be calculated. Otherwise, only p values are calculated using given empirical T stats.
}
}
\details{
This function is to analyze miRNA read count data in two steps: normalization and two-group differential expression test. 

More than one normalization method can be specified in \code{method} when \code{ncpu = 1}. Method \code{GP-Theta} is suggested. There are also other choices of the normalization methods. More details about \code{GP-Quantile}, \code{GP-Theta} can be found in our article.

When sample size is large, \code{maxIter} must be specified (500 by default). Smaller \code{maxIter} may save run time but to be too small may make the empirical distribution unreliable. 

Besides of specifying appropriately small value of \code{maxIter}, it is suggested to make \code{ncpu} larger than 1, where parallel computation is applied. In parallel computing process, permutation task is splitted into parts of almost equal size, each of which will be processed by a core. And the calculation of p values is also paralleled by dividing the miRs into subsets, each of which is processed by a core.
}
\value{
A \code{GPSmle} object. See \link{deGPS_mRNA}.
}
\references{
deGPS: a Powerful and Flexible Framework for Detecting Differential Expression in RNA-Sequencing Studies
}
\seealso{
\link{deGPS_mRNA}
}
\examples{
\dontrun{
##Generate Random samples from GP(theta, lambda)
examData <- newExampleData(nRNA = 100, groupSize = 2, lambda = 0.9, 
theta = 3, ptol = 1e-15)
str(examData)

##Differential Expression Tests for miRNA
examRes <- GPSmle(data = examData$data, group = examData$group, method = "GP-Theta", 
type = "pvalue", ncpu = 1, geneid = paste("G", 1:100, sep = ""))
str(examRes)

topTags(examRes, n = 10, method = "BH")

plot(examRes)
}
}
\keyword{ normalization }
\keyword{ de-test }
