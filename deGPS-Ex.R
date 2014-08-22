pkgname <- "deGPS"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('deGPS')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("GPSmle.default")
### * GPSmle.default

flush(stderr()); flush(stdout())

### Name: GPSmle.default
### Title: Generalized Poisson Statistical Maximum Likelihood Estimation
###   (default)
### Aliases: GPSmle.default
### Keywords: normalization de-test

### ** Examples

## Not run: 
##D ##Generate Random samples from GP(theta, lambda)
##D examData <- newExampleData(nRNA = 100, groupSize = 2, lambda = 0.9, 
##D theta = 3, ptol = 1e-15)
##D str(examData)
##D 
##D ##Differential Expression Tests for miRNA
##D examRes <- GPSmle(data = examData$data, group = examData$group, method = "GP-Theta", 
##D type = "pvalue", ncpu = 1, geneid = paste("G", 1:100, sep = ""))
##D str(examRes)
##D 
##D topTags(examRes, n = 10, method = "BH")
##D 
##D plot(examRes)
## End(Not run)



cleanEx()
nameEx("deGPS-package")
### * deGPS-package

flush(stderr()); flush(stdout())

### Name: deGPS-package
### Title: Normalization and Two-group Differential Expression Test for
###   RNA-seq Data
### Aliases: deGPS-package deGPS
### Keywords: GP RNA-seq DifferentialExpression

### ** Examples

## Not run: 
##D ### Analyze real RNA-seq data
##D ### Compare "Early Embryo" to "Late Embryo" in fly data
##D data(flyData)
##D 
##D i <- 1
##D 
##D group <- rep(1:2, each = 6)
##D simuData <- as.matrix(flyData$data[ , -1])
##D simuData <- simuData[ , flyData$compIdx[[i]]]
##D 
##D 
##D ###remove genes of all-zero read counts
##D simuData <- simuData[apply(simuData, 1, function(x) !all(x == 0)), ]    
##D 
##D ### apply deGPS on fly data with the empirical T stats downloaded at 
##D ### https://www.dropbox.com/s/if5ido5vd8rzff5/empirical.T.stats.fly.RData
##D ### you can set empirical.T.stats as NULL to get your own empirical values
##D ### note that for non-parallelized computation, it may take hours.
##D 
##D load("empirical.T.stats.fly.RData")
##D empirical.T.stats.fly <- empTAll[i]
##D names(empirical.T.stats.fly) <- "GP-Theta"
##D 
##D ### Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
##D flyRes <- deGPS_mRNA(data = simuData, group = group, 
##D ncore = 6, nSubcore = 6, method = "GP-Theta", 
##D empirical.T.stats = empirical.T.stats.fly)
##D 
##D pvalue <- as.vector(flyRes$pvalue)
##D adj.p <- p.adjust(pvalue, method = "BH")
##D topTags(pvalue)
##D 
##D ### Generate Random samples from GP(theta, lambda)
##D examData <- newExampleData(nRNA = 100, groupSize = 6, lambda = 0.9, 
##D theta = 3, ptol = 1e-15)
##D str(examData)
##D 
##D ### Differential Expression Tests
##D examRes <- deGPS_mRNA(data = examData$data, group = examData$group, 
##D method = "GP-Theta", nSubcore = 2, ncore = 2, geneid = paste("G", 1:100, sep = ""))
##D str(examRes)
##D 
##D ### Generate simulated RNA-seq data from compcodeR package
##D require(compcodeR)
##D 
##D samples.per.cond <- 5
##D random.outlier.high.prob <- 0.1
##D n.vars <- 10000
##D 
##D examData <- generateSyntheticData(dataset = "simuData",
##D n.vars = n.vars, samples.per.cond = samples.per.cond, n.diffexp = floor(n.vars * 0.1),
##D repl.id = 1, seqdepth = 1e+07, fraction.upregulated = 0.5,
##D between.group.diffdisp = FALSE, filter.threshold.total = 1,
##D filter.threshold.mediancpm = 0, fraction.non.overdispersed = 0,
##D random.outlier.high.prob = random.outlier.high.prob,
##D output.file = "simuData_repl1.rds")
##D 
##D group <- examData@sample.annotations$condition
##D 
##D ### Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
##D examRes <- deGPS_mRNA(data = examData@count.matrix, group = group, 
##D method = "GP-Theta", nSubcore = 6, ncore = 6, geneid = paste("G", 1:nrow(examData@count.matrix), sep = ""))
##D str(examRes)
##D 
## End(Not run)



cleanEx()
nameEx("deGPS_mRNA")
### * deGPS_mRNA

flush(stderr()); flush(stdout())

### Name: deGPS_mRNA
### Title: Normalization and Two-group Differential Expression Test for
###   mRNA Read Count Data
### Aliases: deGPS_mRNA
### Keywords: deGPS mRNA

### ** Examples

## Not run: 
##D ### See the example in "flyData" for real data analysis and the comparison between deGPS 
##D ### and other widely-used methods
##D 
##D ##Generate Random samples from GP(theta, lambda)
##D examData <- newExampleData(nRNA = 100, groupSize = 6, lambda = 0.9, 
##D theta = 3, ptol = 1e-15)
##D str(examData)
##D 
##D ##Differential Expression Tests
##D examRes <- deGPS_mRNA(data = examData$data, group = examData$group, 
##D method = "GP-Theta", nSubcore = 2, ncore = 2, geneid = paste("G", 1:100, sep = ""))
##D str(examRes)
##D topTags(examRes, n = 10, method = "BH")
##D 
##D ###Generate simulated RNA-seq data from compcodeR package
##D require(compcodeR)
##D 
##D samples.per.cond <- 5
##D random.outlier.high.prob <- 0.1
##D n.vars <- 10000
##D 
##D examData <- generateSyntheticData(dataset = "simuData",
##D n.vars = n.vars, samples.per.cond = samples.per.cond, n.diffexp = floor(n.vars * 0.1),
##D repl.id = 1, seqdepth = 1e+07, fraction.upregulated = 0.5,
##D between.group.diffdisp = FALSE, filter.threshold.total = 1,
##D filter.threshold.mediancpm = 0, fraction.non.overdispersed = 0,
##D random.outlier.high.prob = random.outlier.high.prob,
##D output.file = "simuData_repl1.rds")
##D 
##D group <- examData@sample.annotations$condition
##D 
##D ###Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
##D examRes <- deGPS_mRNA(data = examData@count.matrix, group = group, 
##D method = "GP-Theta", nSubcore = 6, ncore = 6, geneid = paste("G", 1:nrow(examData@count.matrix), sep = ""))
##D str(examRes)
##D topTags(examRes, n = 10, method = "BH")
##D 
## End(Not run)



cleanEx()
nameEx("flyData")
### * flyData

flush(stderr()); flush(stdout())

### Name: flyData
### Title: A real RNA-seq data set of fly
### Aliases: flyData
### Keywords: fly

### ** Examples

## Not run: 
##D ### load required packages
##D 
##D require(edgeR)
##D require(DESeq)
##D require(DESeq2)
##D require(ggplot2)
##D require(gridExtra)
##D require(deGPS)
##D 
##D data(flyData)
##D str(flyData)
##D 
##D ################################################################
##D #### The list of flyData contains:
##D #### data: read counts table with the first row as gene names
##D #### groupInfo: the state name of each sample in the data
##D #### compIdx: six subgroups of indices of samples for analysis
##D ################################################################
##D #### choose the i-th subgroup as an example:
##D #### i = 1: Early vs Late Embryo
##D #### i = 2: Late Embryo vs Larval
##D #### i = 3: Larval vs Adult
##D #### i = 4: Early Embryo vs Larval
##D #### i = 5: Early Embryo vs Adult
##D #### i = 6: Late Embryo vs Adult
##D ################################################################
##D 
##D i <- 1
##D 
##D group <- rep(1:2, each = 6)
##D simuData <- as.matrix(flyData$data[ , -1])
##D simuData <- simuData[ , flyData$compIdx[[i]]]
##D titleName <- names(flyData$compIdx)[i]             ### the name used in the title of the final plot
##D 
##D ###remove genes of all-zero read counts
##D simuData <- simuData[apply(simuData, 1, function(x) !all(x == 0)), ]    
##D 
##D ### apply deGPS on fly data with the empirical T stats downloaded at 
##D ### https://www.dropbox.com/s/if5ido5vd8rzff5/empirical.T.stats.fly.RData
##D ### you can set empirical.T.stats as NULL to get your own empirical values
##D ### note that for non-parallelized computation, it may take hours.
##D 
##D load("empirical.T.stats.fly.RData")
##D empirical.T.stats.fly <- empTAll[i]
##D names(empirical.T.stats.fly) <- "GP-Theta"
##D 
##D ### Make sure you have 6 cores before running deGPS_mRNA with nSubcore = ncore = 6
##D flyRes <- deGPS_mRNA(data = simuData, group = group, 
##D ncore = 6, nSubcore = 6, method = "GP-Theta", 
##D empirical.T.stats = empirical.T.stats.fly)
##D 
##D pvalue <- as.vector(flyRes$pvalue)
##D 
##D ### compare result with edgeR and DESeq, DESeq2
##D 
##D d0 <- DGEList(counts = simuData, group = group)
##D design <- model.matrix(~ group, data = d0$samples)
##D d <- try(calcNormFactors(d0, method = "TMM"))
##D d <- try(estimateGLMCommonDisp(d, design, verbose = TRUE))
##D d <- try(estimateGLMTrendedDisp(d, design))
##D efit <- try(glmQLFTest(d, design, coef = 2))
##D edge2Res <- efit$table$PValue
##D 
##D d <- try(estimateGLMTagwiseDisp(d, design))
##D efit <- try(glmFit(d, design))
##D efit1 <- try(glmLRT(efit, coef = 2))
##D edge1Res <- efit1$table$PValue
##D 
##D cds1 <- newCountDataSet(simuData, group)
##D cds2 <- estimateSizeFactors(cds1)
##D cds3 <- try(estimateDispersions(cds2))
##D if("try-error" ##D 
##D if("try-error" ##D 
##D fitType = "local"))
##D res <- nbinomTest(cds3, 1, 2)
##D deseqRes <- res$pval
##D deseqRes[is.na(deseqRes)] <- 1
##D 
##D dds <- DESeqDataSetFromMatrix(countData = simuData,
##D colData = data.frame(group = group),
##D design = ~ group)
##D dds <- DESeq(dds)
##D res <- results(dds)
##D deseq2Res <- res$pvalue
##D 
##D ### plot the overlap of DEs
##D pAll <- cbind(pvalue, edge1Res, edge2Res, deseqRes, deseq2Res)
##D 
##D pAll[is.na(pAll)] <- 1
##D 
##D pAll <- apply(pAll, 2, p.adjust, method = "BH")
##D 
##D overLapTemp <- overLapTemp1 <- matrix(NA, ncol(pAll), ncol(pAll))
##D 
##D for(ii in 1:ncol(pAll)){
##D for(jj in 1:ncol(pAll)){
##D overLapTemp[ii, jj] <- sum(pAll[ , ii] < 0.05 & pAll[ , jj] < 0.05) / sum(pAll[ , ii] < 0.05)
##D }
##D }
##D 
##D for(ii in 1:ncol(pAll)){
##D for(jj in 1:ncol(pAll)){
##D overLapTemp1[ii, jj] <- sum(pAll[ , ii] < 0.05 & pAll[, jj] < 0.05)
##D }
##D }
##D 
##D compName <- c("deGPS", "edgeR1", "edgeR2", "DESeq", "DESeq2")
##D dimnames(overLapTemp) <- dimnames(overLapTemp1) <- list(compName, compName)
##D 
##D jpeg("flyOverLap.jpg", width = 1600, height = 700)
##D 
##D a <- levelplot(overLapTemp, xlab = "", ylab = "", main = list(paste("Overlap Proportion", titleName), cex = 2), 
##D scale = list(cex = 1.3), 
##D colorkey = list(labels = list(cex = 1.2)), 
##D panel=function(...) {
##D arg <- list(...)
##D panel.levelplot(...)
##D panel.text(rep(1:nrow(overLapTemp), ncol(overLapTemp)), 
##D rep(1:ncol(overLapTemp), each = nrow(overLapTemp)), round(as.vector(overLapTemp), 2), cex = 1.5)}
##D )
##D 
##D b <- levelplot(overLapTemp1, xlab = "", ylab = "", main = list(paste("Overlap Number", titleName), cex = 2), 
##D scale = list(cex = 1.3), 
##D colorkey = list(labels = list(cex = 1.2)), 
##D panel=function(...) {
##D arg <- list(...)
##D panel.levelplot(...)
##D panel.text(rep(1:nrow(overLapTemp1), ncol(overLapTemp1)), 
##D rep(1:ncol(overLapTemp1), each = nrow(overLapTemp1)), as.vector(overLapTemp1), cex = 1.5)}
##D )
##D grid.arrange(a, b, ncol=2)
##D dev.off()
## End(Not run)



cleanEx()
nameEx("newExampleData")
### * newExampleData

flush(stderr()); flush(stdout())

### Name: newExampleData
### Title: Generate example data for GPSmle and deGPS_mRNA
### Aliases: newExampleData
### Keywords: random gp

### ** Examples

## Not run: 
##D ####Different Lambda and Theta for Two Groups
##D examData <- newExampleData(nRNA = 100, groupSize = 2, lambda = c(0.5, 0.9), 
##D theta = c(3, 10), ptol = 1e-15)
##D 
##D ####Same Lambda and Theta for Two Groups
##D examData <- newExampleData(nRNA = 100, groupSize = 2, lambda = 0.9, theta = 3, 
##D ptol = 1e-15)
##D 
##D 
##D ###Generate simulated RNA-seq data from compcodeR package
##D require(compcodeR)
##D 
##D samples.per.cond <- 5
##D random.outlier.high.prob <- 0.1
##D n.vars <- 10000
##D 
##D examData <- generateSyntheticData(dataset = "simuData",
##D n.vars = n.vars, samples.per.cond = samples.per.cond, n.diffexp = floor(n.vars * 0.1),
##D repl.id = 1, seqdepth = 1e+07, fraction.upregulated = 0.5,
##D between.group.diffdisp = FALSE, filter.threshold.total = 1,
##D filter.threshold.mediancpm = 0, fraction.non.overdispersed = 0,
##D random.outlier.high.prob = random.outlier.high.prob,
##D output.file = "simuData_repl1.rds")
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
