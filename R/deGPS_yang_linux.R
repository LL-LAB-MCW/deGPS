#function definition

##################################################################################

toptags <- function(x, n = 10, method = "BH", significance = 0.05){
  if(!("GPSmle" %in% class(x))) x <- list(pvalue = x)
  x$pvalue <- as.matrix(x$pvalue)
  if(ncol(x$pvalue) > 1) stop("multiple rows found in pvalues.")
  adj.p <- p.adjust(x$pvalue, method = method)
  
  if(all(adj.p >= significance)) stop("No DEs in the result.")
  
  pvalueIdx <- adj.p < significance
  if(is.null(dimnames(x$pvalue)[[1]])) geneName <- "no given gene names"
  else geneName <- dimnames(x$pvalue)[[1]][pvalueIdx]
  
  pvalue <- x$pvalue[pvalueIdx]
  
  adj.p <- adj.p[pvalueIdx]
  
  if(sum(pvalueIdx) < n) warning(paste("Only", sum(pvalueIdx), "significant genes found at significant level", significance, "\n"))
  resOrder <- order(adj.p)[1:min(sum(pvalueIdx), n)]
  
  adj.p <- adj.p[resOrder]
  
  pvalue <- x$pvalue[pvalueIdx][resOrder]
  
  geneid <- seq_len(nrow(x$pvalue))[pvalueIdx][resOrder]
  
  resList <- list(pvalue = pvalue, adj.pvalue = adj.p, geneName = geneName, geneID = geneid, significance = significance)
  return(resList)
}

##################################################################################

#Density Function of General Poisson Distribution

dgpois <- function(x, lambda = 0, theta = 1, log = F)
{
  prob <- log(theta) + (x - 1) * log(theta + lambda * x) - theta - lambda * x - lgamma(x + 1)
  if(!log)prob <- exp(prob)
  return(prob)
}

##################################################################################

#Distribution Function of General Poisson Distribution

pgpois <- function(x, lambda = 0, theta = 1, lower.tail = T, log.p = F, tol = 1e6)
{
  tempGP <- fitGP_mRNA(x, tol = 1e-15)
  lambda <- tempGP$lambda
  theta <- tempGP$theta
  
  prob <- cumsum(dgpois(0:min(max(x),tol), lambda = lambda, theta = theta, log = F))
  if(!lower.tail)prob <- 1 - prob
  if(log.p)prob <- log(prob)
  
  x[x>tol] <- tol
  return(prob[x+1])
}

##################################################################################

gpCdf <- function(lambda = 0, theta = 1, maxP = maxP){
  for(i in 1:100){
    if(maxP < 1 - dgpois(10 ^ i, lambda = lambda, theta = theta, log = F)) break
  }
  y <- 0:10^i
  prob <- dgpois(y, lambda = lambda, theta = theta, log = F)
  prob <- cumsum(prob)
  return(prob)
}

##################################################################################

rzigp <- function(n, phi = 0, lambda = 0, theta = 1, ptol = 1e-10){
  signVec <- rbinom(n, size = 1, prob = 1 - phi)
  inflatedZero <- sum(signVec == 0)
  quantileVec <- runif(n - inflatedZero)
  prob <- gpCdf(lambda = lambda, theta = theta, maxP = 1 - ptol)
  ranGen <- unlist(lapply(quantileVec, function(x) which.min(order(c(x, prob))) - 1))
  ranGen1 <- rep(0, n)
  ranGen1[signVec == 1] <- ranGen
  return(ranGen1)
}

##################################################################################

newExampleData <- function(nRNA = 100, groupSize = 5, lambda = 0, theta = 1, ptol = 1e-10){
  warning("The results are just random samples from specified GP distribution, which is FAR AWAY from real RNA-seq data.
          You can use R package 'compcodeR' to generate simulated RNA-seq data. See the examples in the manual for R codes and details.")
  if(any(lambda <= 0)) stop("lambda should be larger than 0.")
  if(any(theta <= 0)) stop("theta should be larger than 0.")
  if(length(lambda) > 2 || length(theta) > 2) stop("lambda or theta cannot have length larger than 2.")
  if(length(lambda) == 1) lambda <- rep(lambda, 2)
  if(length(theta) == 1) theta <- rep(theta, 2)
  genData <- NULL
  for(i in 1:groupSize){
    tempData <- cbind(rzigp(nRNA, phi = 0, lambda = lambda[1], theta = theta[1], ptol = ptol), rzigp(nRNA, phi = 0, lambda = lambda[2], theta = theta[2], ptol = ptol))
    genData <- cbind(genData, tempData)
  }
  group <- rep(1:2, groupSize)
  genData <- genData[ , order(group)]
  return(list(data = genData, group = rep(1:2, each = groupSize)))
}

##################################################################################

preprocess2 <- function (x, data.type = "MAS5", threshold = 1, LOWESS = FALSE, percent = 50) 
{
  x <- as.matrix(na.exclude(x))
  if (data.type == "MAS4" || data.type == "MAS5") {
    x <- quartile.normalize(x, percent = percent)
  }
  if (data.type == "MAS4" || data.type == "MAS5" || data.type == 
        "dChip") {
    if (length(x[x < threshold]) != 0) {
      x[x < threshold] <- threshold
    }
  }
  x <- logb(x, 2)
  if (LOWESS) {
    y <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    y[, 1] <- x[, 1]
    for (i in 2:ncol(x)) {
      y[, i] <- lowess.normalize(x[, 1], x[, i])
    }
    x <- y
  }
  return(x)
}

##################################################################################
### calculate empirical T-test statistic

fun3_unpaired <- function(x, group){
  x1 <- x[group == 1]
  x2 <- x[group == 2]
  if(length(unique(x1)) == 1 && length(unique(x2)) == 1){
    if(mean(x1) != mean(x2))return(1e5)
    else return(0)
  }
  if(length(unique(x1)) == 1)return(abs(mean(x1) - mean(x2))/sd(x2)/sqrt(2 / length(x1)))
  if(length(unique(x2)) == 1)return(abs(mean(x1) - mean(x2))/sd(x1)/sqrt(2 / length(x1)))
  tteststat <- try(t.test(x1, x2)$statistic)
  if("try-error" %in% class(tteststat))return(NA)
  else return(abs(tteststat))
}

##################################################################################
#####normalization
normalFun <- function(dataSim, method = c("Lowess", "GP", "Quantile", "TMM", "GP2")){
  dataSimLowess <- dataSimMLE <- dataSimGP <- dataSimQuantile <- dataSimTMM <- dataSimGP2 <- NULL
  
  if("Lowess" %in% method){
    dataSimLowess <- try(preprocess2(dataSim, LOWESS = T, percent = 90))
    if("try-error" %in% class(dataSimLowess))dataSimLowess <- NULL
  }
  
  if("GP" %in% method){
    dataSimGP <- apply(dataSim, 2, pgpois)
  }
  
  if("Quantile" %in% method){
    dataSimQuantile <- try(normalizeBetweenArrays(dataSim, method = "quantile"))
    if("try-error" %in% class(dataSimQuantile))dataSimQuantile <- NULL
  }
  
  if("TMM" %in% method){
    factorTMM <- try(calcNormFactors(dataSim, method = "TMM"))
    if("try-error" %in% class(factorTMM)){
      dataSimTMM <- NULL
    }else{
      dataSimTMM <- dataSim * rep(factorTMM, each = nrow(dataSim))
    }
  }
  
  if("GP2" %in% method){
    parGP2 <- apply(dataSim, 2, function(x) fitGP_mRNA(x, tol = 1e-15)$theta)
    dataSimGP2 <- dataSim / rep(parGP2, each = nrow(dataSim))
  }
  
  dataNormal <- list(raw = dataSim)
  if(!is.null(dataSimGP))dataNormal <- c(dataNormal, list(GP = dataSimGP))
  if(!is.null(dataSimGP2))dataNormal <- c(dataNormal, list(GP2 = dataSimGP2))
  if(!is.null(dataSimLowess))dataNormal <- c(dataNormal, list(Lowess = dataSimLowess))
  if(!is.null(dataSimQuantile))dataNormal <- c(dataNormal, list(Quantile = dataSimQuantile))
  if(!is.null(dataSimTMM))dataNormal <- c(dataNormal, list(TMM = dataSimTMM))
  return(dataNormal)
}


##################################################################################
#######################Calculate pvalue
calPvalueFun <- function(StatAllArray, dataNormal, 
                         method = c("Lowess", "GP", "Quantile", "TMM", "GP2"), 
                         group = rep(1:2, each = 5), paired = FALSE){
  
  paired <- FALSE
  fun3 <- fun3_unpaired
  
  dataSim <- dataNormal[["raw"]]
  dataSimLowess <- dataNormal[["Lowess"]]
  dataSimGP <- dataNormal[["GP"]]
  dataSimQuantile <- dataNormal[["Quantile"]]
  dataSimTMM <- dataNormal[["TMM"]]
  dataSimGP2 <- dataNormal[["GP2"]]
  
  pValueAll1 <- NULL
  
  if("Lowess" %in% method){
    if(!is.null(dataSimLowess) && !is.null(StatAllArray$Lowess)){
      StatLowess <- apply(dataSimLowess, 1, function(x) fun3(x, group))
      RefCdf <- StatAllArray$Lowess[StatAllArray$Lowess < 1e5]  
      pvalueLowess <- unlist(lapply(StatLowess, function(x)mean(RefCdf > x, na.rm = TRUE)))
    }else{
      pvalueLowess <- rep(NA, nrow(dataSim))
    }
    pValueAll1 <- cbind(pValueAll1, Lowess = pvalueLowess)
  }
  
  if("GP" %in% method){
    StatGP <- apply(dataSimGP, 1, function(x) fun3(x, group))
    RefCdf <- StatAllArray$GP[StatAllArray$GP < 1e5]
    pvalueGP <- unlist(lapply(StatGP, function(x)mean(RefCdf > x, na.rm = TRUE)))
    pValueAll1 <- cbind(pValueAll1, GP = pvalueGP)
  }
  
  if("Quantile" %in% method){
    if(!is.null(dataSimQuantile) && !is.null(StatAllArray$Quantile)){
      StatQuantile <- apply(dataSimQuantile, 1, function(x) fun3(x, group))
      RefCdf <- StatAllArray$Quantile[StatAllArray$Quantile < 1e5]
      pvalueQuantile <- unlist(lapply(StatQuantile, function(x)mean(RefCdf > x, na.rm = TRUE)))
    }else{
      pvalueQuantile <- rep(NA, nrow(dataSim))
    }
    pValueAll1 <- cbind(pValueAll1, Quantile = pvalueQuantile)
  }
  
  if("TMM" %in% method){
    if(!is.null(dataSimTMM) && !is.null(StatAllArray$TMM)){
      StatTMM <- apply(dataSimTMM, 1, function(x) fun3(x, group))
      RefCdf <- StatAllArray$TMM[StatAllArray$TMM < 1e5]
      pvalueTMM <- unlist(lapply(StatTMM, function(x)mean(RefCdf > x, na.rm = TRUE)))
    }else{
      pvalueTMM <- rep(NA, nrow(dataSim))
    }
    pValueAll1 <- cbind(pValueAll1, TMM = pvalueTMM)
  }
  
  if("GP2" %in% method){
    StatGP2 <- apply(dataSimGP2, 1, function(x) fun3(x, group))
    RefCdf <- StatAllArray$GP2[StatAllArray$GP2 < 1e5]  
    pvalueGP2 <- unlist(lapply(StatGP2, function(x)mean(RefCdf > x, na.rm = TRUE)))
    pValueAll1 <- cbind(pValueAll1, GP2 = pvalueGP2)
  }
  pValueAll1 <- as.matrix(pValueAll1)
  #str(pValueAll1)
  pValueAll1 <- pValueAll1[ , method, drop = FALSE]
  return(pValueAll1)
}

#######################Revised Calculate pvalue for large volumn of RNAs

funRevisedP <- function(idx, xRange, ecdf, x){
  rangeIdx <- ((xRange[idx] - 1) * 0.01 * length(ecdf) + 1):(xRange[idx] * 0.01 * length(ecdf))
  mean(ecdf[rangeIdx] > x[idx], na.rm = TRUE) * 0.01 + (100 - xRange[idx]) * 0.01
  #1 - which.min(order(x[idx], ecdf[rangeIdx])) / length(ecdf) + (xRange[idx] - 1) * 0.1
}

calPvalueFunRevised <- function(StatAllArray, dataNormal, 
                                method = c("Lowess", "GP", "Quantile", "TMM", "GP2"), group = rep(1:2, each = 5), 
                                paired = FALSE){
  
  paired <- FALSE
  fun3 <- fun3_unpaired
  
  dataSim <- dataNormal[["raw"]]
  dataSimLowess <- dataNormal[["Lowess"]]
  dataSimGP <- dataNormal[["GP"]]
  dataSimQuantile <- dataNormal[["Quantile"]]
  dataSimTMM <- dataNormal[["TMM"]]
  dataSimGP2 <- dataNormal[["GP2"]]
  
  pValueAll1 <- NULL
  
  if("Lowess" %in% method){
    if(!is.null(dataSimLowess) && !is.null(StatAllArray$Lowess)){
      StatLowess <- apply(dataSimLowess, 1, function(x) fun3(x, group))
      RefCdf <- StatAllArray$Lowess[StatAllArray$Lowess < 1e5] 
      
      RefCdf <- RefCdf[order(RefCdf)]
      quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
      
      xRange <- NULL
      for(range_i in 1:length(StatLowess)){
        tempRange <- which.min(order(c(StatLowess[range_i], quan100)))
        xRange <- c(xRange, tempRange)
      }
      pvalueLowess <- unlist(lapply(1:length(StatLowess), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatLowess))
      
    }else{
      pvalueLowess <- rep(NA, nrow(dataSim))
    }
    pValueAll1 <- cbind(pValueAll1, Lowess = pvalueLowess)
  }
  
  if("GP" %in% method){
    StatGP <- apply(dataSimGP, 1, function(x) fun3(x, group))
    RefCdf <- StatAllArray$GP[StatAllArray$GP < 1e5]
    
    RefCdf <- RefCdf[order(RefCdf)]
    quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
    
    xRange <- NULL
    for(range_i in 1:length(StatGP)){
      tempRange <- which.min(order(c(StatGP[range_i], quan100)))
      xRange <- c(xRange, tempRange)
    }
    pvalueGP <- unlist(lapply(1:length(StatGP), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatGP))
    
    pValueAll1 <- cbind(pValueAll1, GP = pvalueGP)
  }
  
  if("Quantile" %in% method){
    if(!is.null(dataSimQuantile) && !is.null(StatAllArray$Quantile)){
      StatQuantile <- apply(dataSimQuantile, 1, function(x) fun3(x, group))
      RefCdf <- StatAllArray$Quantile[StatAllArray$Quantile < 1e5]
      
      RefCdf <- RefCdf[order(RefCdf)]
      quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
      
      xRange <- NULL
      for(range_i in 1:length(StatQuantile)){
        tempRange <- which.min(order(c(StatQuantile[range_i], quan100)))
        xRange <- c(xRange, tempRange)
      }
      pvalueQuantile <- unlist(lapply(1:length(StatQuantile), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatQuantile))
      
    }else{
      pvalueQuantile <- rep(NA, nrow(dataSim))
    }
    pValueAll1 <- cbind(pValueAll1, Quantile = pvalueQuantile)
  }
  
  if("TMM" %in% method){
    if(!is.null(dataSimTMM) && !is.null(StatAllArray$TMM)){
      StatTMM <- apply(dataSimTMM, 1, function(x) fun3(x, group))
      RefCdf <- StatAllArray$TMM[StatAllArray$TMM < 1e5]
      
      RefCdf <- RefCdf[order(RefCdf)]
      quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
      
      xRange <- NULL
      for(range_i in 1:length(StatTMM)){
        tempRange <- which.min(order(c(StatTMM[range_i], quan100)))
        xRange <- c(xRange, tempRange)
      }
      pvalueTMM <- unlist(lapply(1:length(StatTMM), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatTMM))
      
    }else{
      pvalueTMM <- rep(NA, nrow(dataSim))
    }
    pValueAll1 <- cbind(pValueAll1, TMM = pvalueTMM)
  }
  
  if("GP2" %in% method){
    StatGP2 <- apply(dataSimGP2, 1, function(x) fun3(x, group))
    RefCdf <- StatAllArray$GP2[StatAllArray$GP2 < 1e5]  
    
    RefCdf <- RefCdf[order(RefCdf)]
    quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
    
    xRange <- NULL
    for(range_i in 1:length(StatGP2)){
      tempRange <- which.min(order(c(StatGP2[range_i], quan100)))
      xRange <- c(xRange, tempRange)
    }
    pvalueGP2 <- unlist(lapply(1:length(StatGP2), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatGP2))
    
    pValueAll1 <- cbind(pValueAll1, GP2 = pvalueGP2)
  }
  pValueAll1 <- as.matrix(pValueAll1)
  #str(pValueAll1)
  pValueAll1 <- pValueAll1[ , method, drop = FALSE]
  return(pValueAll1)
}

##################################################################################

##################################################################################
###########################Calculate Ecdf
listCat <- function(total, part){
  listRes <- list()
  length(listRes) <- length(total)
  for(i in 1:length(total)){
    listRes[[i]] <- c(total[[i]], part[[i]])
  }
  return(listRes)
}

ecdfFun <- function(dataNormal, lengthComb, method = "GP2", group, paired = FALSE, combAll, combAll2){
  
  paired <- FALSE
  
  #################################
  ########Specify stat-fun and resampling times for paired and unpaired test
  
  fun3 <- fun3_unpaired
  
  #################################
  ########Get normalized data
  
  dataSim <- dataNormal[["raw"]]
  
  StatAllArray <- list()
  length(StatAllArray) <- length(method)
  
  #################################
  ########Generate Empirical T-stst
  cat("Start randomly shuffling the samples...  \n")
  for(statarray_i in 1:lengthComb)
  {
    ########Generate permutated index
    shufflingIdx <- c(combAll[ , statarray_i], combAll2[ , ncol(combAll) - statarray_i + 1])
    StatAllArrayTemp <- lapply(dataNormal[-1], function(y) apply(y[ , shufflingIdx], 1, function(x) fun3(x, group)))
    StatAllArray <- listCat(StatAllArray, StatAllArrayTemp)
  }
  names(StatAllArray) <- names(dataNormal)[-1]
  return(StatAllArray)
}

ecdfFunMulti <- function(x, maxIter, dataNormal, method = "GP2", group, paired, lengthCombTotal, combAllTotal, combAllTotal2){
  iterIdx <- ((x - 1) * maxIter + 1):min(x * maxIter, lengthCombTotal)
  #iterIdx1 <- max(ncol(combAllTotal2) - lengthCombTotal + 1, ncol(combAllTotal2) - x * maxIter + 1):(ncol(combAllTotal2) - (x - 1) * maxIter)
  combAll <- combAllTotal[ , iterIdx]
  combAll2 <- combAllTotal2[ , ncol(combAllTotal2):1][ , rev(iterIdx)]
  lengthComb <- length(iterIdx)
  a <- ecdfFun(dataNormal = dataNormal, method = method, lengthComb = lengthComb, group = group, paired = paired, combAll = combAll, combAll2 = combAll2)
  return(a)
}


##################################################################################
nameFun <- function(x, method){
  names(x) <- c("raw", method)
  return(x)
}

GPSmleEst <- function(data, group = rep(1:2, each = 5), type = c("normalization", "ecdf", "pvalue"), dataNormal = NULL, 
                      empirical.T.stats = NULL, method = c("Lowess", "GP", "Quantile", "TMM", "GP2"), maxIter = 500, paired = FALSE,
                      ncpu = 1, geneid = NULL){
  i = numeric(1)
  
  method[method == "GP-Quantile"] <- "GP"
  method[method == "GP-Theta"] <- "GP2"
  
  if(is.null(dataNormal)){
    dataNormal <- normalFun(data, method = method)
  }else{
    names(dataNormal)[names(dataNormal) %in% "GP-Quantile"] <- "GP"
    names(dataNormal)[names(dataNormal) %in% "GP-Theta"] <- "GP2"
  }
  
  if(type != "normalization"){
    if(is.null(empirical.T.stats)){
      warning("No given empirical T-stats. They will be calculated by randomly shuffling the samples. \n")
      nGroup1 <- sum(group == 1)
      nGroup2 <- sum(group == 2)
      minNgroup <- round(min(nGroup1, nGroup2) / 2)
      
      combAll <- try(combn(length(group), nGroup1))
      if("try-error" %in% class(combAll)) warning("Too many samples to transverse all the possible permutations. maxIter times of permutations will be applied with default value of 500 for miR.")
      
      theoryMax <- factorial(length(group)) / factorial(nGroup1) / factorial(nGroup2) / (as.integer(nGroup1 == nGroup2) + 1)
      
      if(maxIter > theoryMax){
        maxIter <- theoryMax
        warning("The specified maxIter is larger than the number of all possible permutations of samples. It is re-assigned with that maximum number.")
      }
      if("try-error" %in% class(combAll) || maxIter < theoryMax){ 
        
        #maxIter <- floor(maxIter / ncpu) * ncpu
        lengthComb <- maxIter
        #maxIter <- lengthComb
        
        combAll <- combAll2 <- NULL
        for(comb_i in 1:lengthComb){
          combTemp <- sample(1:length(group), length(group))
          combAll <- cbind(combAll, combTemp[group == 1])
          combAll2 <- cbind(combAll2, combTemp[group == 2])
        }
        combAll2 <- combAll2[ , lengthComb:1]
        
      }else{
        warning("Your specification implies the transversing of all possible permutations.")
        
        if(nGroup1 == nGroup2){
          combAll2 <- combAll
          lengthComb <- ncol(combAll) / 2
        }else{
          combAll2 <- combn(length(group), nGroup2)
          lengthComb <- ncol(combAll)
        }
        
      }
      
      #cat("right multi: ", c(combAll[ , 1], combAll2[ , ncol(combAll)]), "\n")
      
      combAllTotal <- combAll
      combAllTotal2 <- combAll2
      lengthCombTotal <- lengthComb
      
      if(ncpu == 1){
        StatAllArray <- ecdfFun(dataNormal = dataNormal, method = method, lengthComb = lengthComb, group = group, paired = paired, 
                                combAll = combAll, combAll2 = combAll2)
      }else{
        
        if(length(method) != 1) stop("only one normalization method can be used when ncpu is larger than 1.")
        
        cat("Calculating empirical values of T-test stats by randomly shuffling the samples.... \n")
        
        cl <- makeCluster(ncpu)
        registerDoParallel(cl)
        
        # foreach 1
          StatAllArrayMLETemp <- foreach(i = 1:ncpu, .combine = c, .export = c("ecdfFun", "ecdfFunMulti", "fun3_unpaired", "listCat")) %dopar% 
            ecdfFunMulti(x=i, maxIter = ceiling(maxIter / ncpu), dataNormal = dataNormal, method = method, group = group, paired = paired, combAllTotal = combAllTotal,
                         combAllTotal2 = combAllTotal2, lengthCombTotal = lengthCombTotal)        
        stopCluster(cl)
        
        StatAllArray <- list()
        length(StatAllArray) <- 1
        names(StatAllArray) <- method
        StatAllArray[[1]] <- unlist(StatAllArrayMLETemp)
      }
    }else{
      StatAllArray <- empirical.T.stats
      names(StatAllArray)[names(StatAllArray) == "GP-Theta"] <- "GP2"
      names(StatAllArray)[names(StatAllArray) == "GP-Quantile"] <- "GP"
    }
    if(type != "ecdf"){
      if(ncpu > 1){
        
        beginIdx <- seq(1, by = floor(nrow(data) / ncpu), length.out = ncpu)
        endIdx <- c(seq(floor(nrow(data) / ncpu), by = floor(nrow(data) / ncpu), length.out = ncpu - 1), nrow(data))
        dataNormalTemp <- lapply(seq_len(ncpu), function(x) list(dataNormal[[1]][beginIdx[x]:endIdx[x], ], dataNormal[[2]][beginIdx[x]:endIdx[x], ]))
        dataNormalTemp <- lapply(dataNormalTemp, nameFun, method = method)
        
        cat("Start calculating p values...  \n")
        
        cl <- makeCluster(ncpu)
        registerDoParallel(cl)
        
       
          if(length(StatAllArray[[1]]) > 1e6){
		pValueAll1 <- foreach(
				i = 1:ncpu, 
				.combine = c, .export = c("calPvalueFunRevised", "fun3_unpaired", "funRevisedP")
		) %dopar% calPvalueFunRevised(
		            dataNormal =  dataNormalTemp[[i]], StatAllArray = StatAllArray, method = method, group = group, paired = paired
			)
		} else {
			pValueAll1 <- foreach(
				i = 1:ncpu, .combine = c, .export = c("calPvalueFun", "fun3_unpaired")) %dopar% calPvalueFun(
           			dataNormal = dataNormalTemp[[i]], StatAllArray = StatAllArray, method = method, group = group, paired = paired
			)
		}
       
        pValueAll1 <- as.matrix(pValueAll1)
        dimnames(pValueAll1) <- list(NULL, method)
        stopCluster(cl)
        
      }else{
        if(length(StatAllArray[[1]]) > 1e6) pValueAll1  <- calPvalueFunRevised(StatAllArray, dataNormal,  method = method, group = group, paired = paired)
        else  pValueAll1  <- calPvalueFun(StatAllArray, dataNormal,  method = method, group = group, paired = paired)
      }
      
      pValueAll1 <- as.matrix(pValueAll1)
      if(!is.null(geneid)) dimnames(pValueAll1)[[1]] <- geneid
      if(!is.null(dimnames(data)[[1]])) dimnames(pValueAll1)[[1]] <- dimnames(data)[[1]]
      
      names(StatAllArray)[names(StatAllArray) %in% "GP"] <- "GP-Quantile"
      names(StatAllArray)[names(StatAllArray) %in% "GP2"] <- "GP-Theta"
      
      dimnames(pValueAll1)[[2]][dimnames(pValueAll1)[[2]] %in% "GP"] <- "GP-Quantile"
      dimnames(pValueAll1)[[2]][dimnames(pValueAll1)[[2]] %in% "GP2"] <- "GP-Theta"
      
    }else{
      pValueAll1 <- NULL
    }
    
    names(StatAllArray)[names(StatAllArray) %in% "GP"] <- "GP-Quantile"
    names(StatAllArray)[names(StatAllArray) %in% "GP2"] <- "GP-Theta"
    
  }else{
    StatAllArray <- NULL
    pValueAll1 <- NULL
  }
  
  log2FC <- apply(data / rep(apply(data, 2, mean, na.rm = TRUE), each = nrow(data)), 1, function(x) log2(mean(x[group == 2], na.rm = TRUE) / mean(x[group == 1], na.rm = TRUE)))
  log2FC[sum(data[ , group == 1]) == 0] <- Inf
  log2FC[sum(data[ , group == 1]) == 0 && sum(data[ , group == 2]) == 0] <- 0
  
  method[method == "GP"] <- "GP-Quantile"
  method[method == "GP2"] <- "GP-Theta"
  
  names(dataNormal)[names(dataNormal) %in% "GP"] <- "GP-Quantile"
  names(dataNormal)[names(dataNormal) %in% "GP2"] <- "GP-Theta"
  
  resList <- list(normalized.data = dataNormal, log2FoldChange = log2FC, empirical.T.stats = StatAllArray, pvalue = pValueAll1, paired = paired,
                  method = method, type = type)
  return(resList)
}

##################################################################################

GPSmle <- function(data, group = rep(1:2, each = 5), type = c("pvalue", "normalization", "ecdf"), 
                   method = c("GP-Theta", "Lowess", "GP-Quantile", "Quantile", "TMM"), maxIter = 500, paired = FALSE, ncpu = 1, 
                   geneid = NULL, empirical.T.stats = NULL){
  UseMethod("GPSmle")
}

##################################################################################

GPSmle.default <- function(data, group = rep(1:2, each = 5), type = c("pvalue", "normalization", "ecdf"), method = c("GP-Theta", "Lowess", 
                                                                                                                     "GP-Quantile", "Quantile", "TMM"), maxIter = 500, paired = FALSE, ncpu = 1, geneid = NULL, empirical.T.stats = NULL){
  
  method <- match.arg(method, c("GP-Theta", "Lowess", "GP-Quantile", "Quantile", "TMM"))
  type <- match.arg(type, c("pvalue", "normalization", "ecdf"))
  
  #if("TMM" %in% method) require(edgeR)
  #if("Lowess" %in% method) require(LPE)
  #if("Quantile" %in% method) require(limma)
  
  if(missing(group)) stop("group must be specified.")
  if(length(group) != ncol(data))stop("Length of group must equal to ncol of data.")
  if(length(table(group)) != 2) stop("There are more than two groups in the data.")
  tableGroup <- unique(group)
  group1idx <- group == tableGroup[1]
  group2idx <- group == tableGroup[2]
  group[group1idx] <- 1
  group[group2idx] <- 2
  group <- as.integer(group)
  if(tableGroup[1] != "1") warning(paste("rename group ", tableGroup[1], " as 1", sep = ""))
  if(tableGroup[2] != "2") warning(paste("rename group ", tableGroup[2], " as 2", sep = ""))
  
  if(ncpu > 1 && length(group) < 10) {
    warning("ncpu > 1 is ignored because of small sample size. Only one core will be used.")
    ncpu <- 1
  }
  
  if(missing(ncpu))ncpu <- 1
  
  if(!is.null(geneid)){
    if(length(geneid) != nrow(data)) stop("gene-id is with incorrect length.")
    tableGene <- table(geneid)
    if(any(tableGene != 1)) stop("There are duplicates in gene-id. For mRNA data, 
                                 deGPS only process gene-level read counts, i.e., one read count for each gene. If there are replicates, 
                                 use as.matrix to make each replicate a column.")
  }
  
  estRes <- GPSmleEst(data = data, group = group, type = type, method = method, maxIter = maxIter, paired = paired, 
                      ncpu = ncpu, geneid = geneid, empirical.T.stats = empirical.T.stats)
  class(estRes) <- "GPSmle"
  summary(estRes)
  
  return(estRes)
  }

##################################################################################

plot.GPSmle <- function(x, ...){
  if(!("GPSmle" %in% class(x))) stop("Input object must be GPSmle class.")
  if(is.null(x$pvalue)){
    stop("No pvalues to plot.")
  }else{
    warning("Only pvalues plots are generated.")
    if(ncol(x$pvalue) > 1) par(mfrow = c(ceiling(ncol(x$pvalue) / 2), 2))
    else par(mfrow = c(ncol(x$pvalue), 1))
    for(i in 1:ncol(x$pvalue)){
      hist(x$pvalue[, i], breaks = 100, freq = FALSE, main = paste("P-values Distribution of", dimnames(x$pvalue)[[2]][i]), xlab = "Pvalues")
    }
  }
}

summary.GPSmle <- function(object, ...){
  x <- object
  if(!("GPSmle" %in% class(x))) stop("Input object must be GPSmle class.")
  cat("\n\n\n")
  cat("********************************* \n")
  cat("Normalization method: ", x$method, "\n")
  cat("Results contain Normalized Data")
  if(x$type == "ecdf")cat(", Empirical T-test statistics.\n")
  if(x$type == "pvalue" || x$type == "mRNA")cat(", Empirical T-test statistics, P-values.\n")
  
  cat("Tests are ")
  if(!x$paired) cat("not ")
  cat("paired.\n")
  if(x$type == "pvalue" || x$type == "mRNA"){
    cat("\n********************************* \n")
    cat("After adjusting p values by method 'BH', significant DEs are found: \n")
    for(i in 1:ncol(x$pvalue)){
      cat(x$method[i], ": ")
      sigP <- p.adjust(x$pvalue[ , i], method = "BH") < 0.05
      cat(sum(sigP, na.rm = TRUE), "\t(", round(sum(sigP, na.rm = TRUE) / length(sigP) * 100, 2), "%)", "at 0.05 significant level.\n")
    }
  }
  cat("\n********************************* \n")
  cat("\n\n")
}








##################################################################################
##################################################################################
##################################################################################
################################mRNA
fitGP_mRNA <- function(x, tol = 1e-15)
{
  x <- x[!is.na(x)]
  lambda <- uniroot(function(lambda){sum(x * (x - 1) / (mean(x) + (x - mean(x)) * lambda)) - sum(x)}, c(0, 1 - tol), tol = 1e-15)$root
  theta <- mean(x) * (1 - lambda)
  return(list(lambda = lambda, theta = theta))
}


#######################Parallel Function For ECDF

parFun1Other <- function(x, dataNormal = dataNormal, dataSim = dataSim, fun3 = fun3_unpaired, group = group, combAll1 = combAll1, combAll2 = combAll2){
  lengthComb <- ncol(combAll1)
  StatMLE2 <- NULL
  dataSimMLE <- dataNormal[[x$sampleID + 1]][x$beginID:x$endID, ]
  for(statarray_i in 1:lengthComb)
  {
    shufflingIdx <- c(combAll1[ , statarray_i], combAll2[ , ncol(combAll1) - statarray_i + 1]) 
    if(!is.null(dataSimMLE))StatMLE2 <- c(StatMLE2, apply(dataSimMLE[ , shufflingIdx], 1, function(x) fun3(x, group)))
  }
  return(list(res = list(idx = x, ecdf = StatMLE2)))
}

#######################Parallel Function For Pvalue

parFun2Other <- function(x, dataNormal = dataNormal, ecdf = ecdf, dataSim = dataSim, fun3 = fun3_unpaired, group = group){
  
  RefCdf <- ecdf[[x$sampleID]]
  RefCdf <- RefCdf[RefCdf < 1e5]
  RefCdf <- RefCdf[order(RefCdf)]
  quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
  
  dataSimMLE <- dataNormal[[x$sampleID + 1]][x$beginID:x$endID, ]
  StatMLE2 <- apply(dataSimMLE, 1, function(x) fun3(x, group))
  
  xRange <- NULL
  for(range_i in 1:length(StatMLE2)){
    tempRange <- which.min(order(c(StatMLE2[range_i], quan100)))
    xRange <- c(xRange, tempRange)
  }
  if(length(RefCdf) > 1e6) pvalueMLE2 <- unlist(lapply(1:length(StatMLE2), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatMLE2))
  else pvalueMLE2 <- unlist(lapply(StatMLE2, function(x) mean(RefCdf > x, na.rm = TRUE)))
  return(list(res = list(idx = x, pvalue = pvalueMLE2)))
}

#######################Parallel Function For Pvalue regularized t stats
parFun2OtherSamT <- function(x, dataNormal = dataNormal, ecdf = ecdf, dataSim = dataSim, fun3 = fun3_unpaired, group = group){
  
  RefCdf <- ecdf[[x$sampleID]]
  RefCdf <- RefCdf[RefCdf < 1e5]
  RefCdf <- RefCdf[order(RefCdf)]
  quan100 <- quantile(RefCdf, p = seq(0.01, 0.99, 0.01))
  
  dataSimMLE <- dataNormal[[x$sampleID + 1]][x$beginID:x$endID, ]
  #StatMLE2 <- apply(dataSimMLE, 1, function(x) fun3(x, group))
  StatMLE2 <- sam.stat(t(dataSimMLE), group)
  
  xRange <- NULL
  for(range_i in 1:length(StatMLE2)){
    tempRange <- which.min(order(c(StatMLE2[range_i], quan100)))
    xRange <- c(xRange, tempRange)
  }
  if(length(RefCdf) > 1e6) pvalueMLE2 <- unlist(lapply(1:length(StatMLE2), funRevisedP, xRange = xRange, ecdf = RefCdf, x = StatMLE2))
  else pvalueMLE2 <- unlist(lapply(StatMLE2, function(x) mean(RefCdf > x, na.rm = TRUE)))
  return(list(res = list(idx = x, pvalue = pvalueMLE2)))
}

#######################All references for MLE2 normalization
normalFun_mRNA <- function(dataSim, method = c("Lowess", "GP", "Quantile", "TMM", "GP2")){
  dataSimLowess <- dataSimMLE <- dataSimGP <- dataSimQuantile <- dataSimTMM <- dataSimGP2 <- NULL
  method[method == "GP-Quantile"] <- "GP"
  method[method == "GP-Theta"] <- "GP2"
  if("Lowess" %in% method){
    dataSimLowess <- try(preprocess2(dataSim, LOWESS = T, percent = 90))
    if("try-error" %in% class(dataSimLowess))dataSimLowess <- NULL
  }
  
  if("GP" %in% method){
    dataSimGP <- apply(dataSim, 2, pgpois)
  }
  
  if("Quantile" %in% method){
    dataSimQuantile <- try(normalizeBetweenArrays(dataSim, method = "quantile"))
    if("try-error" %in% class(dataSimQuantile))dataSimQuantile <- NULL
  }
  
  if("TMM" %in% method){
    factorTMM <- try(calcNormFactors(dataSim, method = "TMM"))
    if("try-error" %in% class(factorTMM)){
      dataSimTMM <- NULL
    }else{
      dataSimTMM <- dataSim * rep(factorTMM, each = nrow(dataSim))
    }
  }
  
  if("GP2" %in% method){
    parGP2 <- apply(dataSim, 2, function(x) fitGP_mRNA(x)$theta)
    dataSimGP2 <- dataSim / rep(parGP2, each = nrow(dataSim))
  }
  
  dataNormal <- list(raw = dataSim)
  if(!is.null(dataSimGP))dataNormal <- c(dataNormal, list(GP = dataSimGP))
  if(!is.null(dataSimGP2))dataNormal <- c(dataNormal, list(GP2 = dataSimGP2))
  if(!is.null(dataSimLowess))dataNormal <- c(dataNormal, list(Lowess = dataSimLowess))
  if(!is.null(dataSimQuantile))dataNormal <- c(dataNormal, list(Quantile = dataSimQuantile))
  if(!is.null(dataSimTMM))dataNormal <- c(dataNormal, list(TMM = dataSimTMM))
  return(dataNormal)
}

######################Parallel Main Function with Given nCore and nSubcore
####################################Other Methods

mRnaEcdfOther <- function(dataSim = dataSim, dataNormal = dataNormal, ecdf = NULL, group = rep(1:2, each = 5), method = "GP2", 
                          nSubcore = 4, ncore = 4, paired = FALSE, maxIter = 150, 
                          mixFDR = FALSE){
  
  method[method == "GP-Quantile"] <- "GP"
  method[method == "GP-Theta"] <- "GP2"
  
  paired <- FALSE
  fun3 <- fun3_unpaired
  
  sampleID <- rep(1:length(method), each = nSubcore)
  eachITER <- floor(nrow(dataSim) / nSubcore)
  if(eachITER < 1){
    warning("nSubcore is excess, only part of it will be used.")
    nSubcore <- nrow(dataSim)
    eachITER <- 1
  }
  
  beginID <- rep(seq(1, by = eachITER, length.out = nSubcore), length(method))
  endID <- rep(c(seq(eachITER, by = eachITER, length.out = nSubcore - 1), nrow(dataSim)), length(method))
  listIdx <- lapply(1:(length(method) * nSubcore), function(x) list(sampleID = sampleID[x], beginID = beginID[x], endID = endID[x]))
  #str(listIdx)
  resAll <- NULL
  
  nGroup1 <- sum(group == 1)
  nGroup2 <- sum(group == 2)
  
  if(is.null(ecdf)){
    
    theoryMax <- factorial(length(group)) / factorial(nGroup1) / factorial(nGroup2) / (as.integer(nGroup1 == nGroup2) + 1)
    
    if(maxIter > theoryMax){
      maxIter <- theoryMax
      warning("The specified maxIter is larger than the number of all possible permutations of samples. It is re-assigned with that maximum number.")
    }
    
    if(length(group) > 10){
      
      combAll1 <- combAll2 <- NULL
      for(statarray_i in 1:maxIter){
        shufflingIdx <- sample(1:length(group), length(group))
        combAll1 <- cbind(combAll1, shufflingIdx[group == 1])
        combAll2 <- cbind(combAll2, shufflingIdx[group == 2])
      }
      combAll2 <- combAll2[ , ncol(combAll2):1]
      
    }else{
      warning("Your specification implies the transversing of all possible permutations.")
      if(nGroup1 == nGroup2){
        combAll <- combn(1:(2 * nGroup1), nGroup1)
        combAll1 <- combAll[ , 1:(ncol(combAll) / 2)]
        combAll2 <- combAll[ , (ncol(combAll) / 2 + 1):ncol(combAll)]
      }else{
        combAll1 <- combn(1:length(group), nGroup1)
        combAll2 <- combn(1:length(group), nGroup2)
      }
      
    }
    
    cat("Calculating empirical values of T-test stats by randomly shuffling the samples.... \n")
    
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    
     # foreach 3
      resAll <- foreach(i = 1:length(listIdx), .combine = c, .export = c("parFun1Other")) %dopar% 
        parFun1Other(listIdx[[i]], dataNormal = dataNormal, dataSim = dataSim, combAll1 = combAll1, combAll2 = combAll2, fun3 = fun3, group = group)
    
    stopCluster(cl)
    
    ecdfAll <- list()
    length(ecdfAll) <- length(method)
    names(ecdfAll) <- method
    #str(resAll)
    for(i in 1:length(resAll)){
      resTemp <- resAll[[i]]
      #str(resTemp)
      ecdfAll[[resTemp$idx$sampleID]] <- c(ecdfAll[[resTemp$idx$sampleID]], resTemp$ecdf)
    }
  }else{
    names(ecdf)[names(ecdf) == "GP-Theta"] <- "GP2"
    names(ecdf)[names(ecdf) == "GP-Quantile"] <- "GP"
    ecdfAll <- ecdf
  }
  
  ##################Pvalue
  cat("Start calculating p values...  \n")
  if(mixFDR){
    resAll <- lapply(ecdfAll, function(x) fdrtool(c(apply(dataNormal[[2]], 1, fun3_unpaired, group = group), x), 
                                                  statistic = "normal", plot = FALSE))
    pvalueAll <- matrix(1, nrow(dataSim), length(method))
    for(i in 1:length(resAll)){
      pvalueAll[ , i] <- resAll[[i]]$pval[1:nrow(dataSim)]
    }
  }else{
    resAll <- NULL
    
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    
      # foreach 4	
      resAll <- foreach(i = 1:length(listIdx), .combine = c, .export = c("parFun2Other", "funRevisedP")) %dopar% parFun2Other(listIdx[[i]], dataNormal = dataNormal, ecdf = ecdfAll, dataSim = dataSim, fun3 = fun3, group = group)
    
    stopCluster(cl)
  
    pvalueAll <- matrix(1, nrow(dataSim), length(method))
    for(i in 1:length(resAll)){
      resTemp <- resAll[[i]]
      pvalueAll[resTemp$idx$beginID:resTemp$idx$endID, resTemp$idx$sampleID] <- resTemp$pvalue
    }
  }
  pvalueAll <- as.matrix(pvalueAll)
  dimnames(pvalueAll) <- list(NULL, NULL)
  dimnames(pvalueAll)[[2]] <- names(ecdfAll) <- names(dataNormal)[-1]
  return(list(pvalue = pvalueAll, StatAllArray = ecdfAll))
}

####### for regularized t stats

mRnaEcdfOtherSamT <- function(dataSim = dataSim, dataNormal = dataNormal, ecdf = NULL, group = rep(1:2, each = 5), method = "GP2", 
                              nSubcore = 4, ncore = 4, paired = FALSE, maxIter = 150, 
                              mixFDR = FALSE){
  
  method[method == "GP-Quantile"] <- "GP"
  method[method == "GP-Theta"] <- "GP2"
  
  paired <- FALSE
  fun3 <- fun3_unpaired
  
  sampleID <- rep(1:length(method), each = nSubcore)
  eachITER <- floor(nrow(dataSim) / nSubcore)
  if(eachITER < 1){
    warning("nSubcore is excess, only part of it will be used.")
    nSubcore <- nrow(dataSim)
    eachITER <- 1
  }
  
  beginID <- rep(seq(1, by = eachITER, length.out = nSubcore), length(method))
  endID <- rep(c(seq(eachITER, by = eachITER, length.out = nSubcore - 1), nrow(dataSim)), length(method))
  listIdx <- lapply(1:(length(method) * nSubcore), function(x) list(sampleID = sampleID[x], beginID = beginID[x], endID = endID[x]))
  #str(listIdx)
  resAll <- NULL
  
  nGroup1 <- sum(group == 1)
  nGroup2 <- sum(group == 2)
  
  if(is.null(ecdf)){
    
    theoryMax <- factorial(length(group)) / factorial(nGroup1) / factorial(nGroup2) / (as.integer(nGroup1 == nGroup2) + 1)
    
    if(maxIter > theoryMax){
      maxIter <- theoryMax
      warning("The specified maxIter is larger than the number of all possible permutations of samples. It is re-assigned with that maximum number.")
    }
    
    if(length(group) > 10){
      
      combAll1 <- combAll2 <- NULL
      for(statarray_i in 1:maxIter){
        shufflingIdx <- sample(1:length(group), length(group))
        combAll1 <- cbind(combAll1, shufflingIdx[group == 1])
        combAll2 <- cbind(combAll2, shufflingIdx[group == 2])
      }
      combAll2 <- combAll2[ , ncol(combAll2):1]
      
    }else{
      warning("Your specification implies the transversing of all possible permutations.")
      if(nGroup1 == nGroup2){
        combAll <- combn(1:(2 * nGroup1), nGroup1)
        combAll1 <- combAll[ , 1:(ncol(combAll) / 2)]
        combAll2 <- combAll[ , (ncol(combAll) / 2 + 1):ncol(combAll)]
      }else{
        combAll1 <- combn(1:length(group), nGroup1)
        combAll2 <- combn(1:length(group), nGroup2)
      }
      
    }
    
    cat("Calculating empirical values of T-test stats by randomly shuffling the samples.... \n")
    
    #cl <- makeCluster(ncore)
    #registerDoParallel(cl)
    
    #if(.Platform$OS.type == "windows"){
    #resAll <- foreach(i = 1:length(listIdx), .combine = c, .export = c("parFun1Other")) %dopar% 
    #parFun1Other(listIdx[[i]], dataNormal = dataNormal, dataSim = dataSim, combAll1 = combAll1, combAll2 = combAll2, fun3 = fun3, group = group)
    #}else{
    #resAll <- do.call(c, mclapply(listIdx, parFun1Other, mc.cores = ncore, dataNormal = dataNormal, dataSim = dataSim, combAll1 = combAll1, combAll2 = #combAll2, fun3 = fun3, group = group))
    #}
    
    #stopCluster(cl)
    
    ecdfAll <- list()
    length(ecdfAll) <- length(method)
    names(ecdfAll) <- method
    
    #for(i in 1:length(resAll)){
    #resTemp <- resAll[[i]]
    
    #ecdfAll[[resTemp$idx$sampleID]] <- c(ecdfAll[[resTemp$idx$sampleID]], resTemp$ecdf)
    #}
    
    lengthComb <- ncol(combAll1)
    StatMLE2 <- NULL
    dataSimMLE <- dataNormal[[2]]
    for(statarray_i in 1:lengthComb)
    {
      shufflingIdx <- c(combAll1[ , statarray_i], combAll2[ , ncol(combAll1) - statarray_i + 1]) 
      if(!is.null(dataSimMLE))StatMLE2 <- c(StatMLE2, sam.stat(t(dataSimMLE[ , shufflingIdx]), group))
    }
    ecdfAll[[1]] <- StatMLE2
    names(ecdfAll) <- "GP2" 
    
  }else{
    names(ecdf)[names(ecdf) == "GP-Theta"] <- "GP2"
    names(ecdf)[names(ecdf) == "GP-Quantile"] <- "GP"
    ecdfAll <- ecdf
  }
  
  ##################Pvalue
  cat("Start calculating p values...  \n")
  resAll <- NULL
  if(mixFDR){
    resAll <- lapply(ecdfAll, function(x) fdrtool(c(sam.stat(t(dataNormal[[2]]), group), x), 
                                                  statistic = "normal", plot = FALSE))
    pvalueAll <- matrix(1, nrow(dataSim), length(method))
    for(i in 1:length(resAll)){
      pvalueAll[ , i] <- resAll[[i]]$pval[1:nrow(dataSim)]
    }
  }else{
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    
    
      resAll <- foreach(i = 1:length(listIdx), .combine = c, .export = c("parFun2OtherSamT", "funRevisedP")) %dopar% parFun2OtherSamT(listIdx[[i]], dataNormal = dataNormal, ecdf = ecdfAll, dataSim = dataSim, fun3 = fun3, group = group)
    
    stopCluster(cl)
    
    pvalueAll <- matrix(1, nrow(dataSim), length(method))
    for(i in 1:length(resAll)){
      resTemp <- resAll[[i]]
      pvalueAll[resTemp$idx$beginID:resTemp$idx$endID, resTemp$idx$sampleID] <- resTemp$pvalue
    }
  }
  pvalueAll <- as.matrix(pvalueAll)
  dimnames(pvalueAll) <- list(NULL, NULL)
  dimnames(pvalueAll)[[2]] <- names(ecdfAll) <- names(dataNormal)[-1]
  return(list(pvalue = pvalueAll, StatAllArray = ecdfAll))
}


###########################################################################################
###########################################################################################
###########################################################################################
deGPS_mRNA <- function(data, dataNormal = NULL, empirical.T.stats = NULL, group = rep(1:2, each = 5), 
                       method = "GP-Theta", nSubcore = 4, ncore = 4,
                       paired = FALSE, regularized = FALSE, mixFDR = FALSE, maxIter = 150, 
                       geneid = NULL){
  
  method <- match.arg(method, c("GP-Theta", "Lowess", "GP-Quantile", "Quantile", "TMM"))
  
  #if("TMM" %in% method) require(edgeR)
  #if("Lowess" %in% method) require(LPE)
  #if("Quantile" %in% method) require(limma)
  
  paired <- FALSE
  fun3 <- fun3_unpaired
  
  if((regularized | mixFDR) & length(method) > 1) stop("cannot apply regularized T stats on normalization methods other than GP-Theta.")
  if(!regularized) warning("traditional t stats are used when option regularized is FALSE")
  if(!mixFDR) warning("no adjustments made for pooled t stats when calculating p values as 
                      option mixFDR is FALSE.")
  if(missing(group)) stop("group must be specified.")
  if(length(group) != ncol(data))stop("Length of group must equal to ncol of data.")
  if(length(table(group)) != 2) stop("There are more than two groups in the data.")
  tableGroup <- unique(group)
  group1idx <- group == tableGroup[1]
  group2idx <- group == tableGroup[2]
  group[group1idx] <- 1
  group[group2idx] <- 2
  group <- as.integer(group)
  if(tableGroup[1] != "1") warning(paste("rename group ", tableGroup[1], " as 1", sep = ""))
  if(tableGroup[2] != "2") warning(paste("rename group ", tableGroup[2], " as 2", sep = ""))
  
  if(!is.null(geneid)){
    if(length(geneid) != nrow(data)) stop("gene-id is with incorrect length.")
    tableGene <- table(geneid)
    if(any(tableGene != 1)) stop("There are duplicates in gene-id. For mRNA data, 
deGPS only process gene-level read counts, i.e., one read count for each gene. If there are replicates, 
use as.matrix to make each replicate a column.")
  }
  
  if(is.null(dataNormal))dataNormal <- normalFun_mRNA(data, method = method)
  if(regularized) mRnaEcdfRes1 <- mRnaEcdfOtherSamT(dataSim = data, dataNormal = dataNormal, ecdf = empirical.T.stats,
                                                    group = group, method = method, nSubcore = nSubcore, ncore = ncore, paired = paired, maxIter = maxIter, mixFDR = mixFDR)
  else mRnaEcdfRes1 <- mRnaEcdfOther(dataSim = data, dataNormal = dataNormal, ecdf = empirical.T.stats,
                                     group = group, method = method, nSubcore = nSubcore, ncore = ncore, paired = paired, maxIter = maxIter, mixFDR = mixFDR)
  pValueAll1 <- mRnaEcdfRes1$pvalue
  StatAllArray <- mRnaEcdfRes1$StatAllArray
  
  log2FC <- apply(data / rep(apply(data, 2, mean, na.rm = TRUE), each = nrow(data)), 1, function(x) log2(mean(x[group == 2], na.rm = TRUE) / mean(x[group == 1], na.rm = TRUE)))
  log2FC[sum(data[ , group == 1]) == 0] <- Inf
  log2FC[sum(data[ , group == 1]) == 0 && sum(data[ , group == 2]) == 0] <- 0
  
  names(StatAllArray)[names(StatAllArray) == "GP"] <- "GP-Quantile"
  names(StatAllArray)[names(StatAllArray) == "GP2"] <- "GP-Theta"
  
  names(dataNormal)[names(dataNormal) == "GP"] <- "GP-Quantile"
  names(dataNormal)[names(dataNormal) == "GP2"] <- "GP-Theta"
  
  pValueAll1 <- as.matrix(pValueAll1)
  if(!is.null(geneid)) dimnames(pValueAll1)[[1]] <- geneid
  if(!is.null(dimnames(data)[[1]])) dimnames(pValueAll1)[[1]] <- dimnames(data)[[1]]
  
  dimnames(pValueAll1)[[2]][dimnames(pValueAll1)[[2]] == "GP"] <- "GP-Quantile"
  dimnames(pValueAll1)[[2]][dimnames(pValueAll1)[[2]] == "GP2"] <- "GP-Theta"
  
  method[method == "GP"] <- "GP-Quantile"
  method[method == "GP2"] <- "GP-Theta"
  
  resList <- list(normalized.data = dataNormal, log2FoldChange = log2FC, empirical.T.stats = StatAllArray, pvalue = pValueAll1, paired = paired,
                  method = method, type = "mRNA")
  class(resList) <- "GPSmle"
  summary(resList)
  return(resList)
}
