\name{flyData}
\alias{flyData}
\title{
A real RNA-seq data set of fly
}
\description{
A real RNA-seq data set of fly
}
\usage{
data(flyData)
}
\details{
A RNA-seq data set of fly downloaded at \href{http://bowtie-bio.sourceforge.net/recount}{Fly Data}, study \code{modencodefly}
}
\references{
The developmental transcriptome of Drosophila melanogaster, Graveley BR and etc., Nature 2011 Mar 24;471(7339):473-9.
}
\examples{
\dontrun{
### load required packages

require(edgeR)
require(DESeq)
require(DESeq2)
require(ggplot2)
require(gridExtra)
require(deGPS)

data(flyData)
str(flyData)

################################################################
#### The list of flyData contains:
#### data: read counts table with the first row as gene names
#### groupInfo: the state name of each sample in the data
#### compIdx: six subgroups of indices of samples for analysis
################################################################
#### choose the i-th subgroup as an example:
#### i = 1: Early vs Late Embryo
#### i = 2: Late Embryo vs Larval
#### i = 3: Larval vs Adult
#### i = 4: Early Embryo vs Larval
#### i = 5: Early Embryo vs Adult
#### i = 6: Late Embryo vs Adult
################################################################

i <- 1

group <- rep(1:2, each = 6)
simuData <- as.matrix(flyData$data[ , -1])
simuData <- simuData[ , flyData$compIdx[[i]]]
titleName <- names(flyData$compIdx)[i]             ### the name used in the title of the final plot

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

### If you cannot download the RData on dropbox, please run following codes instead
### flyRes <- deGPS_mRNA(data = simuData, group = group, 
### ncore = 6, nSubcore = 6, method = "GP-Theta", 
### empirical.T.stats = NULL)

pvalue <- as.vector(flyRes$pvalue)

### compare result with edgeR and DESeq, DESeq2

d0 <- DGEList(counts = simuData, group = group)
design <- model.matrix(~ group, data = d0$samples)
d <- try(calcNormFactors(d0, method = "TMM"))
d <- try(estimateGLMCommonDisp(d, design, verbose = TRUE))
d <- try(estimateGLMTrendedDisp(d, design))
efit <- try(glmQLFTest(d, design, coef = 2))
edge2Res <- efit$table$PValue

d <- try(estimateGLMTagwiseDisp(d, design))
efit <- try(glmFit(d, design))
efit1 <- try(glmLRT(efit, coef = 2))
edge1Res <- efit1$table$PValue

cds1 <- newCountDataSet(simuData, group)
cds2 <- estimateSizeFactors(cds1)
cds3 <- try(estimateDispersions(cds2))
if("try-error" %in% class(cds3))cds3 <- try(estimateDispersions(cds2, method = "pooled"))
if("try-error" %in% class(cds3))cds3 <- try(estimateDispersions(cds2, method = "pooled", 
fitType = "local"))
res <- nbinomTest(cds3, 1, 2)
deseqRes <- res$pval
deseqRes[is.na(deseqRes)] <- 1

dds <- DESeqDataSetFromMatrix(countData = simuData,
colData = data.frame(group = group),
design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
deseq2Res <- res$pvalue

### plot the overlap of DEs
pAll <- cbind(pvalue, edge1Res, edge2Res, deseqRes, deseq2Res)

pAll[is.na(pAll)] <- 1

pAll <- apply(pAll, 2, p.adjust, method = "BH")

overLapTemp <- overLapTemp1 <- matrix(NA, ncol(pAll), ncol(pAll))

for(ii in 1:ncol(pAll)){
for(jj in 1:ncol(pAll)){
overLapTemp[ii, jj] <- sum(pAll[ , ii] < 0.05 & pAll[ , jj] < 0.05) / sum(pAll[ , ii] < 0.05)
}
}

for(ii in 1:ncol(pAll)){
for(jj in 1:ncol(pAll)){
overLapTemp1[ii, jj] <- sum(pAll[ , ii] < 0.05 & pAll[, jj] < 0.05)
}
}

compName <- c("deGPS", "edgeR1", "edgeR2", "DESeq", "DESeq2")
dimnames(overLapTemp) <- dimnames(overLapTemp1) <- list(compName, compName)

jpeg("flyOverLap.jpg", width = 1600, height = 700)

a <- levelplot(overLapTemp, xlab = "", ylab = "", main = list(paste("Overlap Proportion", titleName), cex = 2), 
scale = list(cex = 1.3), 
colorkey = list(labels = list(cex = 1.2)), 
panel=function(...) {
arg <- list(...)
panel.levelplot(...)
panel.text(rep(1:nrow(overLapTemp), ncol(overLapTemp)), 
rep(1:ncol(overLapTemp), each = nrow(overLapTemp)), round(as.vector(overLapTemp), 2), cex = 1.5)}
)

b <- levelplot(overLapTemp1, xlab = "", ylab = "", main = list(paste("Overlap Number", titleName), cex = 2), 
scale = list(cex = 1.3), 
colorkey = list(labels = list(cex = 1.2)), 
panel=function(...) {
arg <- list(...)
panel.levelplot(...)
panel.text(rep(1:nrow(overLapTemp1), ncol(overLapTemp1)), 
rep(1:ncol(overLapTemp1), each = nrow(overLapTemp1)), as.vector(overLapTemp1), cex = 1.5)}
)
grid.arrange(a, b, ncol=2)
dev.off()
}
}
\keyword{fly}

