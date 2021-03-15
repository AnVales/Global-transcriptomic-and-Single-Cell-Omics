# PRACTICE 2 EXERCISES #
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


# PART 1 #
tablecsv2 <-read.csv("./TableS5.csv",row.names=1) #path
colnames(tablecsv2)
dim(tablecsv2)
class(tablecsv2)

matrix_2<-as.matrix(tablecsv2)
boxplot(matrix_2)

# PART 2 + PART 5#

PCA2_1<-princomp(tablecsv2,cor=FALSE,scores=TRUE)
plot(PCA2_1)

summary(PCA2_1) 
# Comp1: Standard deviation     212.0878090
# Comp2: Standard deviation     148.9216057
# Comp3: Standard deviation     137.2364961

variance = PCA2_1$sdev[1:9]/sum(PCA2_1$sdev[1:9])
scores = PCA2_1$scores
loadings = PCA2_1$loadings

plot(loadings)
text(loadings,labels=row.names(loadings),pos=c(2,2,3,1,1,1,1,1,1,1,1,1,1,1))

df12 = data.frame(PC1=loadings[,1], PC2=loadings[,2])
df34 = data.frame(PC3=loadings[,3], PC4=loadings[,4])
rownames(df12) <-rownames(loadings)
rownames(df34) <-rownames(loadings)

df12s = data.frame(PC1=scores[,1], PC2=scores[,2])
df34s = data.frame(PC3=scores[,3], PC4=scores[,4])
rownames(df12s) <-rownames(scores)
rownames(df34s) <-rownames(scores)

# loadings
plot(df12$PC1, df12$PC2)
text(df12,labels=row.names(df12),pos=c(2,2,3,1,1,1,1,1,1,1,1,1,1,1))

plot(df34$PC3, df34$PC4)
text(df34,labels=row.names(df34),pos=c(2,2,3,1,1,1,1,1,1,1,1,1,1,1))

# scores
plot(df12s)

plot(df34s)

# PART 3 #
# X<-subset(CellTypeData, E30>1&CO2<1&COR<1)
subset_marcadores<-subset(tablecsv2, WOX5>1&APL<1&CO2<1&COBL9<1&COR<1&E30<1&GL2<1&PET111<1&S17<1&S18<1&S32<1&S4<1&SCR<1&WER<1&WOL<1)
subset_marcadores_matrix<-as.matrix(subset_marcadores)
heatmap(subset_marcadores_matrix, col = newcolors)

# PART 4 #
PCA2_2<-princomp(subset_marcadores,cor=FALSE,scores=TRUE)
plot(PCA2_2, main="varianza")

variance2 = PCA2_2$sdev[1:9]/sum(PCA2_2$sdev[1:9])
scores2 = PCA2_2$scores
loadings2 = PCA2_2$loadings

summary(PCA2_2) 

df12_2 = data.frame(PC1=loadings2[,1], PC2=loadings2[,2])
df34_2 = data.frame(PC3=loadings2[,3], PC4=loadings2[,4])
rownames(df12_2) <-rownames(loadings2)
rownames(df34_2) <-rownames(loadings2)

df12_2s = data.frame(PC1=scores2[,1], PC2=scores2[,2])
df34_2s = data.frame(PC3=scores2[,3], PC4=scores2[,4])
rownames(df12_2s) <-rownames(scores2)
rownames(df34_2s) <-rownames(scores2)

# loadings
plot(df12_2,xlim=c(-1.1,0.1), main= "Plot loadings PC1 vs. PC2")
text(df12_2,labels=row.names(df12_2),pos=c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))

plot(df34_2, main = "Plot loadings PC3 vs. PC4")
text(df34_2,labels=row.names(df34_2),pos=c(2,2,3,1,1,1,1,1,1,1,1,1,1,1))

# scores
plot(df12_2s, main= "Scores PC1 vs. PC2")

plot(df34_2s,main = "Scores PC3 vs. PC4")

# Plot variance

plot(variance2, main= "Scores PC1 vs. PC2")

# Comp1: Standard deviation     212.0878090
# Comp2: Standard deviation     148.9216057
# Comp3: Standard deviation     137.2364961



