\name{toptags}
\alias{toptags}
\title{
The top significant genes or miRs
}
\description{
The top significant genes or miRs
}
\usage{
toptags(x, n = 10, method = "BH", significance = 0.05)
}
\arguments{
  \item{x}{
A \code{GPSmle} object returned by \code{\link{GPSmle}} or \code{\link{deGPS_mRNA}} or a p value vector. If it is a \code{GPSmle} object, only one column of p values is allowed in this function.
}
  \item{n}{
the number of required top significant genes or miRs.
}
  \item{method}{
the adjust method of multiple testing p values, same as R function \code{p.adjust}.
}
  \item{significance}{
the significant level of the DE.
}
}
\details{
This function is to find the significant genes or miRs at the given significant level. If you want to call this function, make sure that \code{library(deGPS)} is called after other packages, such as \code{edgeR}, since those packages also contain function named \code{toptags}.
}
\value{
A list containing:
  \item{pvalue}{
the ordered p values of significant genes or miRs. 
}
  \item{adj.pvalue}{
the ordered adjusted p values by specified method. 
}
  \item{method}{
p value adjusted method, same as the method in R function \code{p.adjust}
}
  \item{geneName}{
names of the genes or miRs.
}
  \item{geneid}{
the row indice of the genes or miRs in given p values.
}
  \item{significance}{
the specified significant level.
}
}


