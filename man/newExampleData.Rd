\name{newExampleData}
\alias{newExampleData}
\title{
Generate example data for GPSmle and deGPS_mRNA
}
\description{
Randomly generate Generalized Poisson distributed samples with given theta and lambda.  
}
\usage{
newExampleData(nRNA = 100, groupSize = 5, lambda = 0, theta = 1, ptol = 1e-10)
}
\arguments{
  \item{nRNA}{
The number of genes or miRs.
}
  \item{groupSize}{
A integer represents The number of samples in each group. Note that the function can only generate two equal gropus.
}
  \item{lambda}{
The lambda parameter of GP distribution. The values must be within (0, 1). Since miRNA/mRNA read counts tend to be overdispersed, we constain the lambda of example data larger than zero to be similar to the real cases. Note that the length of \code{lambda} can be either one or two, representing the equal or unequal lambda for each group.
}
  \item{theta}{
The theta parameter of GP distribution. Must be positive values. It can be a single value or a numeric vector with length two.
}
  \item{ptol}{
The tolerance of probabilites of GP distribution. The default value is 1e-15, since regular R can not tell the difference smaller than 1e-15. See \code{details} for more explanations.
}
}
\details{
The resulting data set contains two GP distributed groups with the specified lambda and theta as the parameters, with each group containing \code{groupSize} samples. Note that the length of \code{lambda} and \code{theta} can be either one or two, representing the same or different GP for two groups. Moreover, the data is neither H0 (non-DE) nor H1 (with DE) data. Every element is a random sample of the given GP and one can not tell whether one single row in the data is DE. However, with large amount of RNAs in the data, it tends to be H0 data, since the variability in one particular row is caused by random assignment of two GP distributions. 

The random samples of the specified GP distribution are generated from a multinomial distribution, the domain of which is from zero to a maximum value -- \code{gpMax} + 1. The larger \code{gpMax} is, the closer two distributions are. 

The maximum integer \code{gpMax} is determined as the minimum integer satisfying \eqn{P(x = gpMax) \ge ptol}. The probability for each value from zero to gpMax + 1 is then calculated as the probability of that in specified GP distribution. Note that the probability of \eqn{P(x = gpMax + 1) = 1 - P(x = 0) -\ldots - P(x = gpMax)}. Hence, the smaller \code{ptol} is, the closer the approximated multinomial distribution is to the specified GP distribution. And since regular R, i.e. without particular package which enables more precise calculations, can not tell differences smaller than 1e-15, \code{ptol} is set as 1e-15.

Another way to generate simulated RNA-seq data is the \code{compcodeR} package. Details can be found in the package manual and user guide document. See a simple example in the following example session.
}
\value{
  \item{group}{The group indicator of the resulting data.}
  \item{data}{The resulting data with GP distributon.}
}

\seealso{
\code{\link{deGPS_mRNA}}, \code{\link{GPSmle}}
}
\examples{
\dontrun{
####Different Lambda and Theta for Two Groups
examData <- newExampleData(nRNA = 100, groupSize = 2, lambda = c(0.5, 0.9), 
theta = c(3, 10), ptol = 1e-15)

####Same Lambda and Theta for Two Groups
examData <- newExampleData(nRNA = 100, groupSize = 2, lambda = 0.9, theta = 3, 
ptol = 1e-15)


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
}
}

\keyword{ random }
\keyword{ gp }
