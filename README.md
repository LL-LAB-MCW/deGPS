#deGPS 
|  |            |
| ------------- | ----------- |
|**Version** | 1.0 |
|**Author** | Chen Chu |
|**Depends** | R (>= 3.0.0), **foreach**, doParallel |
|**Suggests** | LPE, limma, edgeR |
|**License** | GPL-2 |
| **Binary source** | <a href="https://github.com/LL-LAB-MCW/deGPS-source-file/blob/master/deGPS_1.0.zip?raw=true"><img src="https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/githublogo.png"></a> <a href="https://degps-rna-seq.googlecode.com/svn/deGPS_1.0.zip"><img src="https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/Google-logo.png"></a>  |
| **Source** | <a href="https://github.com/LL-LAB-MCW/deGPS-source-file/blob/master/deGPS.linux_1.0.tar.gz?raw=true"><img src="https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/githublogo.png"></a> <a href="https://degps-rna-seq.googlecode.com/svn/deGPS.linux_1.0.tar.gz"><img src="https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/Google-logo.png"></a>  |
| **Manual** | <a href="https://github.com/LL-LAB-MCW/deGPS-source-file/blob/master/deGPS-manual.pdf?raw=true"><img src="https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/githublogo.png"></a> <a href="https://degps-rna-seq.googlecode.com/svn/deGPS-manual.rar"><img src="https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/Google-logo.png"></a>  |
| **Following sessions**| **Installation, Details, FAQs, Contact, Reference**|

This package is proposed to analyze RNA-seq data in two steps: normalization and permutation-based differential expression test. 

The permutation step may become time-consuming in large sample size context or in mRNA read count data analysis. Parallel computation is introduced in the main functions to deal with the computational burden. Note that in situations where parallel computation is not necessary, to force parallelization may result in more run time.

Note that deGPS can only handle DE tests between two groups, but the novel GP-based normalization methods can be applied on any RNA-seq data. More Details about the GP-based normalization method and the advantages and limitations of our DE test can be found in our article.

The package also contains function to generate GP distributed data to be an example of RNA-seq data. It has to be pointed out that the simulated data in this package is FAR AWAY from real RNA-seq data, it is just simple GP distributed samples. For more appropriately simulated RNA-seq, **compcodeR** package is suggested, or, alternatively, you can generate H0 and H1 data from real data. Simple examples are given in the manual. More R codes referring to the real data based simulation or **compcodeR** based simulation can be found in the supplementary materials of our article.

Please do not hesitate to contact the author if you have any questions or find any flaws in the package.

##Installation
####deGPS
- **via devtools**
 ```
  install.packages("devtools")
  library(devtools)
  install_github("LL-LAB-MCW/deGPS")           ## Windows
  install_github("LL-LAB-MCW/deGPS.linux")     ## Linux
 ```

- **via source (Windows)**
 ```
  install.packages("yourZipFilePath\\deGPS_1.0.zip")  ## change '\' into '\\' in windows path
 ```

 Or you can click the menu of Rgui:

 ![Alt text](https://raw.githubusercontent.com/LL-LAB-MCW/deGPS-source-file/master/ccpic1.jpg?raw=TRUE)

- **via source (Linux)**

 ```
  R CMD INSTALL -l your_R_library_path_in_linux deGPS.linux_1.0.tar.gz
 ```

####Dependency
 ```
  ### Depends
  install.packages(c("foreach", "doParallel"))    ## Windows
  install.packages(c("parallel", "doParallel"))   ## Linux
  
  ### Suggests
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("edgeR", "LPE", "limma"))            ## you can also add "DESeq" & "DESeq2" for comparison
 ```

##Details
- **Normalization**

 New normalization methods based on generalized Poisson distribution, namely **GP-Theta** and  **GP-Quantile**, are contained in the main analysis functions.

 Other popular normalization methods, such as TMM, LOWESS, Quantile, is also available in the package. 

 **GP-Theta**: recommended. Divide every sample by its fitted theta parameter of GP distribution.

 **GP-Quantile**: map every read count to the probability of the fitted GP distribution.


- **Differential expression test**

 deGPS differential expression test is designed for RNA-seq data of small sample size. And it is proved by simulations to successfully control Type I Errors and FDR than widely-used methods, such as edgeR or DESeq, DESeq2.

 P-values are calculated according to the empirical distributions of T-statistic derived from the permutation of normalized RNA-seq data. Reliable empirical distributions can be obtained by pooling T-stat of all miRs or genes together.

 When sample size is too large to iterate through all possible permutations, one can specify **maxIter** to control the maximum permutation times or **ncpu** / **ncore** / **nSubcore** to apply parallel computation. 


- **Comparison with other methods**

 Comprehensive simulations of comparisons have been conducted and the results are presented in our article. The R codes for real data based and **compcodeR** simulated data based  simulation can be found in the supplementary of our article.

##FAQs
- **Q1**: There are errors occurred when running **deGPS**, why DOES NOT it stop?
 
 **A1**: R function **try** is introduced in the source codes of **deGPS**, e.g.:
 ```
 try(combn(length(group), number.per.group))
 ```
 is used in permutation step to determine whether a random sampling can be applied (instead of iterating through all the possible permutations) when **maxIter** is not specified. "try-error" may occur when group size is large. However, it won't affect the permutation, since in this case randomly shuffling will be applied in permutation and the result of **combn** is not needed.
 
 
- **Q2**: Why do you have two different sources for Windows and Linux?
 
 **A2**: Different packages are needed for parallel computation on these two platforms. **parallel** package on Linux and **foreach** for Windows. In practice, sometimes **foreach** fails to do parallelized jobs, instead it just do normal iterations on just one core. For more stable utilization, we choose the function **mcapply** in **parallel** on Linux.
 
 
- **Q3**: Why do you have two different sources for Windows and Linux?
 
 **A3**: Different packages are needed for parallel computation on the two platforms: **parallel** package on Linux and **foreach** for Windows. In practice, sometimes **foreach** fails to do parallelized jobs, instead it just do normal iterations on just one core. For more stable utilization, we choose the function **mcapply** in **parallel** on Linux.
 
 
- **Q4**: Why hasn't your article been published or even accepted?
 
 **A4**: Apparently this is the one asked by myself, and sorry I don't know either. But I think that's the main reason why I have the time to make up this FAQ. lol.
 
##Contact
 ylu@mcw.edu

 chuchen.blueblues@gmail.com

##Reference
 deGPS: a Powerful and Flexible Framework for Detecting Differential Expression in RNA-Sequencing Studies

