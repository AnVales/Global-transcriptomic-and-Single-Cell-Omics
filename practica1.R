# PRACTICE 1 EXERCISES #
if(!require(matlib)){install.packages("matlib")} 
if(!require(rlang)){install.packages("rlang")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(openxlsx)){install.packages("Seurat")}
if(!require(openxlsx)){install.packages("dplyr")}
library(matlib)
library(rlang)
library(openxlsx)
library(Seurat)
library(Matrix)
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")

BiocManager::install("igraph")


install.packages("remotes")

remotes::install_version("SDMTools", "1.1-221")

install.packages('devtools')
devtools::install_version(package = 'Seurat', version = package_version('3.1.1'))

library(Seurat)

# colors
newcolors<-colorRampPalette(colors=c("gray48","skyblue","yellow"))(256)


# PCA among ground tissue samples for shr mutant
tablecsv <-read.csv("./table.csv",row.names=1) #path
colnames(tablecsv)
dim(tablecsv)
class(tablecsv)

# the first 8 columns
data2 <- tablecsv[,1:8]
rownames(data2) <-rownames(tablecsv)
colnames(data2) <- c("shr","WT","shr+BLJ","shr+JKD","shr+MGP","shr+NUC","shr+IME","shr+SCR")
print(data2)
class(data2)

# PCA
PCA1<-princomp(data2,cor=FALSE,scores=TRUE)
plot(PCA1)

PCA1$loadings 
plot(PCA1$loadings)
text(PCA1$loadings,labels=row.names(PCA1$loadings),pos=c(4,4,4,2,2,4))

PCA1$scores
plot(PCA1$scores)
text(PCA1$scores,labels=row.names(PCA1$scores),pos=c(4,4,4,2,2,4))

# create intermediate transcriptome between shr mutant and the wild type (WT)
# transform -> quarter = 0.25*WT+ 0.75*shr.
data2 <- transform(data2, quarter=(0.25*(WT)+0.75*shr), half=(0.50*(WT)+0.50*shr), threequarters=(0.75*(WT)+0.25*shr))

PCA2<-princomp(data2,cor=FALSE,scores=TRUE)
plot(PCA2)

PCA2$loadings 
plot(PCA2$loadings)
text(PCA2$loadings,labels=row.names(PCA2$loadings),pos=c(4,4,4,2,2,4))

PCA$scores
plot(PCA$scores)


# Add the transcriptome of cells corresponding to SCR domain and recalculate PCA.
# Plot variance for components and loadings.
# What might you conclude for several of these transcription factors?
data3<-cbind(data2, tablecsv$SCRdomain )
colnames(data3) <- c("shr","WT","shr+BLJ","shr+JKD","shr+MGP","shr+NUC","shr+IME","shr+SCR","quarter", "half", "threequarters", "SCRdomain")
rownames(data3) <-rownames(tablecsv)

PCA3<-princomp(data3,cor=FALSE,scores=TRUE)
plot(PCA3)

PCA3$loadings 
plot(PCA3$loadings)
text(PCA3$loadings,labels=row.names(PCA3$loadings),pos=c(4,4,4,2,2,4))

PCA3$scores
plot(PCA3$scores)


# Find the most important genes which contributing to observed transcriptomic changes by extracting genes (~20) with highest and lowest score values
# Investigate expression patterns of these genes across original samples using heatmap ().
# Remember functions sort(), head() and tail() might be useful.
sorted_comp1<-sort(PCA3$scores[,1])
top20_1<-head(sorted_comp1, 20)
less20_1<-tail(sorted_comp1,20)

sorted_comp2<-sort(PCA3$scores[,2])
top20_2<-head(sorted_comp2, 20)
less20_2<-tail(sorted_comp2,20)

# Unimos top_20
top_20_t=c(top20_1,top20_2)
top_20_m=as.matrix(top_20_t)
top20=data3[rownames(top_20_m),]
top_20_fm=as.matrix(top20)
heatmap(top_20_fm, col = newcolors)

# Unimos less_20
less_20_t=c(less20_1,less20_2)
less_20_m=as.matrix(less_20_t)
less20=data3[rownames(less_20_m),]
less_20_fm=as.matrix(less20)
heatmap(less_20_fm, col = newcolors)
