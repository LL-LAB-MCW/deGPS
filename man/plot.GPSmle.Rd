\name{plot.GPSmle}
\alias{plot.GPSmle}
\title{
Plot GPSmle
}
\description{
Plot the histograms of GPSmle results.
}
\usage{
\method{plot}{GPSmle}(x, ...)
}
\arguments{
  \item{x}{
the object returned by \code{GPSmle}.
}
  \item{\dots}{
the parameters of \code{plot}
}
}
\details{
See \code{\link{GPSmle.default}} for more details of the output of \code{GPSmle}.
}
\value{
The output depends on the specification of \code{type} in \code{\link{GPSmle}}. If type = "normalization", the histograms of normalized data sets are returned. So are the "ecdf" and "pvalue" or "mRNA" in \code{\link{deGPS_mRNA}}.
}