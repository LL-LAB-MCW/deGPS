\name{deGPS_mRNA}
\alias{deGPS_mRNA}
\title{
Normalization and Two-group Differential Expression Test for mRNA Read Count Data
}
\description{
Normalization and Two-group Differential Expression Test for mRNA Read Count Data
}
\usage{
deGPS_mRNA(data, dataNormal = NULL, empirical.T.stats = NULL, 
group = rep(1:2, each = 5), method = "GP-Theta", nSubcore = 4, 
ncore = 4, paired = FALSE, regularized = FALSE, mixFDR = FALSE, maxIter = 150, geneid = NULL)
}
\arguments{
  \item{data}{
a matrix containing mRNA gene-level read count data. The column represents samples while the row represents genes. Biologcial or technical replicates must be made as new columns in the data.
}
  \item{dataNormal}{
not used anymore.
}
  \item{empirical.T.stats}{
A list of empirical T statistics with names as normalization methods, the format of which must be the same as the \code{empirical.T.stats} in the returned list of \code{deGPS_mRNA} or \code{GPSmle}, i.e., a list of empirical T-stats with the method names as the list name. If null, the empirical T statistics will be calculated. Otherwise, only p values are calculated using given empirical T stats.
}
  \item{group}{
The group indicator. The length of group must equal to ncol(data).
}
  \item{method}{
the methods of normalization. It can be a single charactor or a charactor vector, the values can be "Lowess", "GP-Quantile", "Quantile", "TMM" or "GP-Theta" .
}
  \item{nSubcore}{
The parallel computation strategy splits rows of the data, i.e. mRNAs, into \code{nSubcore} parts. The empirical T-stats are calculated for each part of the data, and so are the p values. The total number of cores needed in the computation is \code{nSubcore} * \code{nMethod}.
}
  \item{ncore}{
The total cpu cores used for the calculations. The specified \code{ncore} can be less than the total number of cores needed (i.e. \code{nSubcore} * \code{nMethod}), in which case the cores will be used repeatedly. Apparently, the maximum utilization of the \code{ncore} cores is reached when \code{mod(nSubcore * nMethod / ncore) == 0}. You may just make \code{ncore = nSubcore} to ensure the efficiency of the parallel computation.
}
  \item{paired}{
The current version of deGPS only contain unpaired test.
}
  \item{regularized}{
Wether regularized t stats are applied instead of traditional t stats. If TRUE, \code{sam.stat} function in package \code{st} will be applied to get empirical distribution.
}
  \item{mixFDR}{
If TRUE, fdrtool is applied on original pooled T statistics to estimate the null empirical distribution and to get the p values. Note that it may be time-consuming if mixFDR is TRUE, since no parallel computation is applied in the process.
}
  \item{maxIter}{
The default value of maxIter is 150. When sample size is large, instead of transversing every possible permutation, randomly shuffling is applied for \code{maxIter} times to obtain the empirical distributions. Larger \code{maxIter} costs longer run time. Note that \code{maxIter} is forced to be not larger than \code{permutationTimes}.
}
  \item{geneid}{
Gene id of the specified data. Biological or technical replicates must be new columns in the data, i.e., duplicates in \code{geneid} are not allowed.
}
}
\details{
This function is to analysis mRNA gene-level read count data in two steps: nomalization and permutation based differential expression test.

More than one normalization method can be specified in one run, method \code{GP-Theta} is suggested. 

In permutation DE test, p values are calculated according to the empirical T statistics obtained by randomly shuffling the samples. You can also specify your own empirical T-stats (or the one you get from another run of the function, which is useful in real data based simulations) in argument \code{empirical.T.stats}. Note that regularized T stats can be used instead of traditional T stats when \code{regularized} is TRUE. Also \code{fdrtool} function in package \code{fdrtool} can be applied on original pooled T stats to estimate the null empirical distribution when \code{mixFDR} is TRUE.

To deal with the burden of computation, e.g. when SampleSize > 10, parallel computing is embeded in the function. By specify \code{nSubcore} and \code{ncore}, parallel computation is applied in the calculations. The abundant genes are divided into \code{nSubcore} subsets, for each of which the empirical T stats are therefore calculated parallelly. The calculation of p values are also parallelly applied on the subsets of mRNAs using empirical T stats obtained by binding \code{nSubcore} subsets. \code{nSubcore * nMethod} cores are needed in total, where \code{nMethod} represents the number of applied methods. 

If the number of cores needed in parallel computing process is larger than \code{ncore}, cores are iteratively used by the introduced R function \code{foreach} (Windows) or \code{mcapply} (Linux). The maximum utilization of the \code{ncore} cores is reached when \code{mod(nSubcore * nMethod / ncore) == 0}. 

Besides of specifying large \code{nSubcore} and \code{ncore}, specify smaller \code{maxIter} can be also useful to make the function more efficiency. Note that \code{maxIter} should not be too small, where the empirical distribution may not be reliable.
}
\value{
A \code{GPSmle} object is returned.
  \item{normalized.data}{A list of normalized data, each element represents one speicified normalization method.}
  \item{log2FoldChange}{The logrithm of fold change}
  \item{empirical.T.stats}{A list of the empirical T-stats of normalied data of different normalization methods, generated by permutation of samples. The length of the T-stats is \code{nRNA} * \code{min(permuationTimes, maxIter)}. \eqn{permutationTimes = \left(\begin{array}{ll} nSample \\ nSample / 2 \end{array}\right) / 2}, if each group has equal size. \eqn{permutationTimes = \left(\begin{array}{l} nSample \\ groupSize \end{array}\right)}, if each group has unequal size.}
  \item{pvalue}{The resulting pvalues. Note that the pvalues may be slightly different in different runs of deGPS for the same data when not all possible permutations are transversed in the calculations of empirical T stats.
}
  \item{paired}{FALSE}
  \item{method}{the normalization methods applied to get the result.}
  \item{type}{"mRNA"}
}
\author{
Chen Chu
}
\seealso{
\code{\link{GPSmle}}
}
\keyword{ deGPS }
\keyword{ mRNA }
\examples{
\dontrun{
### See the example in "flyData" for real data analysis and the comparison between deGPS 
### and other widely-used methods

##Generate Random samples from GP(theta, lambda)
examData <- newExampleData(nRNA = 100, groupSize = 6, lambda = 0.9, 
theta = 3, ptol = 1e-15)
str(examData)

##Differential Expression Tests
examRes <- deGPS_mRNA(data = examData$data, group = examData$group, 
method = "GP-Theta", nSubcore = 2, ncore = 2, geneid = paste("G", 1:100, sep = ""))
str(examRes)
toptags(examRes, n = 10, method = "BH")

###Generate simulated RNA-seq data from compcodeR package
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

###Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
examRes <- deGPS_mRNA(data = examData@count.matrix, group = group, 
method = "GP-Theta", nSubcore = 6, ncore = 6, geneid = paste("G", 1:nrow(examData@count.matrix), sep = ""))
str(examRes)
toptags(examRes, n = 10, method = "BH")

}
}
