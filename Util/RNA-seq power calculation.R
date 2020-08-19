install.packages("/Software/qvalue_2.2.2.tar.gz",repos = NULL, type="source")
library(qvalue)
install.packages("/Software/ssizeRNA_1.2.4.tar.gz",repos = NULL, type="source")
library(ssizeRNA)
set.seed(2016) 
nGenes = 20000# Total number of genes 
pi0 = 0.8#Proportion of non-DE genes
m = 200#psudo number of samples
mu = 10#Average read count for each gene in control group
disp = 0.1#Dispersion parameter for each gene
logfc = log(2)
fdr = 0.05
power = 0.8#Desired power
maxN = 20 #The maximum sample size used for power calculations
size1 <- ssizeRNA_single(nGenes = 20000, pi0 = 0.8, m = 200, mu = 10,
                                        disp = 0.1, logfc = log(2), fdr = 0.05, power = 0.8, maxN = 20)