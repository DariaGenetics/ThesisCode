#Note - project in 'ThesisCode' repository

### libraries
install.packages(...)
library(data.table)
library(ggplot2)
library(SeqArray)
library("Hmisc")
library(FactoMineR)
library(factoextra)
library("writexl")
library(corrplot)
library("ggpubr")
library(dplyr)
library(tidyr)
library(ggfortify)
library(reshape2) #devtools::install_github("hadley/reshape")
library(missMDA)

### load
setwd("~/Desktop")
pheno <- readRDS("/Users/dgg/Desktop/wideform.phenotypedata (1).RDS")
write_xlsx(pheno,"/Users/dgg/Desktop\\AllPhenoNoMeans.xlsx")

'~~~Take data, subtract columns, calculate column means~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
Phenowide <- readRDS("/Users/dgg/Desktop/wideform.phenotypedata (1).RDS") #201 rows; threshold for column inclusion is 100 (50%) non-NA values in rows
AllPhenoNoRalIdWide <- subset(Phenowide, select = -c(ral_id))
w <- as.vector(colSums(is.na(AllPhenoNoRalIdWide)) < 66) # 50% = 100 rows; 33% = 66 rows
Phenonarrow <- AllPhenoNoRalIdWide[, ..w]
AllPhenonarrow <- Phenonarrow %>%
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))

STwide <-fread("/Users/dgg/Desktop/ST-Lines.csv") #231 columns, 158 rows; threshold for column inclusion is 79 (50%) non-NA values in rows
STwide <- subset(STwide, select = -c(V232))
x <- as.vector(colSums(is.na(STwide)) < 52) #50% = 79 rows, 30% = 52 rows
STnarrow <- STwide[, ..x]
STData <- STnarrow %>%
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) #maintained 100 columns

Mixedwide <- fread("/Users/dgg/Desktop/INVST-Lines.csv") #231 columns, 25 rows; threshold for column inclusion is 13 (~50%) non-NA values in rows
Mixedwide <- subset(Mixedwide, select = -c(V232))
y <- as.vector(colSums(is.na(Mixedwide)) < 13) #50% = 13 rows
STnarrow <- STwide[, ..y]
MixedData <- STnarrow %>%
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) #maintained 213 columns


INVwide <- fread("/Users/dgg/Desktop/INV-Lines.csv") #231 columns, 18 rows; threshold for column inclusion is 6 (33.3%) non-NA values in rows
INVwide <- subset(INVwide, select = -c(V232))
z <- as.vector(colSums(is.na(INVwide)) < 6)
INVnarrow <- INVwide[, ..z]
INVData <- INVnarrow %>%
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) #maintained 79 columns

'~~~Correlations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
AllCor <- cor(AllPhenonarrow)
STCor <- cor(STData)
MixedCor <- cor(MixedData)
INVCor <- cor(INVData) #In cor(INVData) : the standard deviation is zero

'~~~Priminitive Heatmaps~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
heatmap(AllCor) #<- just using the correlation for all data

AllPhenoScaled <- scale(AllPhenonarrow)
AllPhenoScaledCor <- cor(AllPhenoScaled)
Heatmap <- heatmap(AllPhenoScaledCor)

'~~~Modulated Modularity Clustering Heatmaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
library("pheatmap")
library("dendsort")

pheatmap(AllCor,fontsize = 3) #<- using non-scaled data (CORRECT)

pheatmap(AllCor,fontsize = 3, kmeans_k = 30) #<- interesting, provides a customizeable horizontal number of vertices


callback = function(hc,...){dendsort(hc)}

pheatmap(AllCor,fontsize = 3, clustering_callback = callback)
STCorHeatmap <- pheatmap(STCor,fontsize = 2, clustering_callback = callback)
INVCorHeatmap <- pheatmap(INVCor,fontsize = 2, clustering_callback = callback) #<- there are NAs

'~~~Creating Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
STScaled <- scale(STData) #normalizes with a mean of 0 and a variance of 1 (standardization)
distSTScaled <- dist(STScaled)
ClustSTScaled <- hclust(distSTScaled) # forms a hierarchical cluster of the data points based on distance metric on dataset objects
plot(ClustSTScaled, ylab = "Height", xlab="Distance")

MixedScaled <- scale(MixedData) #normalizes with a mean of 0 and a variance of 1 (standardization)
distMixedScaled <- dist(MixedScaled)
ClustMixedScaled <- hclust(distMixedScaled) # forms a hierarchical cluster of the data points based on distance metric on dataset objects
plot(ClustMixedScaled, ylab = "Height", xlab="Distance")

INVScaled <- scale(INVData) #normalizes with a mean of 0 and a variance of 1 (standardization)
distINVScaled <- dist(INVScaled)
ClustINVScaled <- hclust(distINVScaled) # forms a hierarchical cluster of the data points based on distance metric on dataset objects
plot(ClustINVScaled, ylab = "Height", xlab="Distance")

'~~~Modulated Modularity Clustering - Net ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
library("modMax")
library("igraph")

#AllCor - plan - make adjacency matrix from the correlation matrix, then use this adjacency matrix to create a new heatmap
AllPhenoScaledCor[upper.tri(AllPhenoScaledCor)] <- 42 # replacing upper triangle of matrix with number that can be used for filtering

ExperimentAllPhenoScaledCor <- melt(AllPhenoScaledCor)

ExperimentAllPhenoScaledCor2 <- filter(ExperimentAllPhenoScaledCor, value != 42) %>% filter(Var1 != Var2)
my_adj_list <- ExperimentAllPhenoScaledCor2 %>% filter(value > 0.60)
dim(my_adj_list)

AllCorNet <- graph.data.frame(my_adj_list, directed = FALSE)
orig_mar <- par()$mar
par(mar=rep(.1, 4))

AllCorCeb <- cluster_edge_betweenness(AllCorNet) # shows correlations between phenotypes, pretty(ier)
plot(AllCorCeb, AllCorNet)

plot(net, layout = layout_components(net), edge.width = E(net)$weight) #vertexshape = none / shows all points with names, ugly

#https://davetang.org/muse/2017/03/16/matrix-to-adjacency-list-in-r/

##removing repeats (done at filter > 0.95)
Lifespan_DR-DietContained05percentYeast_F
Lifespan_Restricted-DR_F

CentroidSize_standard_F
WingCentroidSize_standard_female

WingCentroidSize_standard_male
ControidSize_standard_M

InternocularDistance_standard_male
InternocularDistance_standard_M

InternocularDistance_standard_female
InternocularDistance_standard_F

Lifespan_AL-standardW-5percentYeastAlteration_F
Lifespan_AL-DietCOntained5percentYeast_F

ResistanceToDDT_2HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F #<- really weird formatting
ResistanceToDDT_3HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F

ResistanceToDDT_8HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F #<- really weird formatting
ResistanceToDDT_9HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F
ResistanceToDDT_10HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F
ResistanceToDDT_11HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F
ResistanceToDDT_12HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F
ResistanceToDDT_13HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F
ResistanceToDDT_14HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F
ResistanceToDDT_15HourExposureTo2001.25Acetone-DDTSolutionConcentration05Î1/4g-ml_F

#INVCor - plan - make adjacency matrix from the correlation matrix, then use this adjacency matrix to create a new heatmap
INVCor[upper.tri(INVCor)] <- 42 # replacing upper triangle of matrix with number that can be used for filtering
AdjINVCor <- melt(INVCor)
AdjINVCor2 <- filter(AdjINVCor, value != 42) %>% filter(Var1 != Var2)
AdjListINVCor <- AdjINVCor2 %>% filter(value > 0.80)

INVCorNet <- graph.data.frame(AdjListINVCor, directed = FALSE)
orig_mar <- par()$mar
par(mar=rep(.1, 4))

INVCorCeb <- cluster_edge_betweenness(INVCorNet) # shows correlations between phenotypes, pretty(ier)
plot(INVCorCeb, INVCorNet)

'~~~PCA tables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
##Creating PCA plots
##    name               description
## 1  "$eig"             "eigenvalues"
## 2  "$var"             "results for the variables"
## 3  "$var$coord"       "coord. for the variables"
## 4  "$var$cor"         "correlations variables - dimensions"
## 5  "$var$cos2"        "cos2 for the variables"
## 6  "$var$contrib"     "contributions of the variables"
## 7  "$ind"             "results for the individuals"
## 8  "$ind$coord"       "coord. for the individuals"
## 9  "$ind$cos2"        "cos2 for the individuals"
## 10 "$ind$contrib"     "contributions of the individuals"
## 11 "$call"            "summary statistics"
## 12 "$call$centre"     "mean of the variables"
## 13 "$call$ecart.type" "standard error of the variables"
## 14 "$call$row.w"      "weights for the individuals"
## 15 "$call$col.w"      "weights for the variables"

library(ggplot2)
#finds the scree plot for each phenotype

pcaST <- prcomp(STData, scale.=TRUE)
var_explained_STpca <- data.frame(pcaST= paste0("PC",1:158),
                               var_explained=(pcaST$sdev)^2/sum((pcaST$sdev)^2))
var_explained_STpca <- var_explained_STpca[1:9,1:2]
var_explained_STpca %>%
   ggplot(aes(x=pcaST,y=var_explained, group=1))+
   geom_point(size=4)+
   geom_line()+
   labs(title="Scree plot: PCA on scaled ST data")


pcaINV <- prcomp(INVData, scale.=TRUE)
var_explained_INVpca <- data.frame(pcaINV= paste0("PC",1:18),
                               var_explained=(pcaINV$sdev)^2/sum((pcaINV$sdev)^2))
var_explained_INVpca <- var_explained_INVpca[1:9,1:2]
var_explained_INVpca %>%
   ggplot(aes(x=pcaINV,y=var_explained, group=1))+
   geom_point(size=4)+
   geom_line()+
   labs(title="Scree plot: PCA on scaled INV data")

#finds phenotypes contributing to specific principle components
STDataPCA = PCA(STData, scale.unit=TRUE, ncp=5, graph=T)
res.ST <- PCA(STData, graph = FALSE)
fviz_contrib(res.ST, choice = "var", axes = 1, top = 10) #~ explains the variation from each phenotype

MixedDataPCA = PCA(MixedData, scale.unit=TRUE, ncp=5, graph=T)
res.Mixed <- PCA(MixedData, graph = FALSE)
fviz_contrib(res.Mixed, choice = "var", axes = 1, top = 10) #~ explains the variation from each phenotype

INVDataPCA = PCA(INVData, scale.unit=TRUE, ncp=5, graph=T)
res.INV <- PCA(INVData, graph = FALSE)
fviz_contrib(res.INV, choice = "var", axes = 1, top = 10) #~ explains the variation from each phenotype

#PCA for the whole dataset
library('shinyaframe')
library(ggplot2)

setwd("~/Desktop")
AllPhenoWInversionStatus <- fread("/Users/dgg/Desktop/AllPheno.csv", fill = TRUE)
AllPhenoWInversionStatus <- subset(AllPhenoWInversionStatus, select = -c(V233, V234, V235))


RalID <- as.vector(unlist(AllPhenoWInversionStatus$ral_id))
In2LtStatus <- as.vector(unlist(AllPhenoWInversionStatus$In2Lt))


AllPhenoWInversionStatus <- subset(AllPhenoWInversionStatus, select = -c(ral_id)) #colnumber [1]
AllPhenoWOInversionStatus = subset(AllPhenoWInversionStatus, select = -c(In2Lt) ) #colnumber [231]

AllPhenoPCA1 = imputePCA(AllPhenoWOInversionStatus, scale.unit=TRUE, ncp=5, graph=T) #note - figure out how to impute means as the column means for the INV and ST groups separately
PrincipleComponents <- AllPhenoPCA1$fittedX
PrincipleComponents2 <- as.data.frame(PrincipleComponents) #note - V1 = principle component 1, V2 = principle component 2, etc...
PrincipleComponents2$RalId <- RalID
PrincipleComponents2$In2LtStatus <- In2LtStatus

'~~~~~~Intercept and Inversion Status Graphs~~~~~~~~~~~~~~~~~~~~~~~~'
#incercept and inversions status flux graphs( i = 1:4, 1:10, 1:20, 1:40)
write_xlsx(PrincipleComponents,"/Users/dgg/Desktop/PrincipleComponents.xlsx")

PrincipleComponentsTrimmed <- fread("/Users/dgg/Desktop/PrincipleComponentsWOLastFourRows.csv") #extra 3 rows removed
PrincipleComponentsTrimmed <- fread("/Users/alanbergland/Documents/GitHub/ThesisCode/PrincipleComponentsWOLastFourRows.csv")
MeltedPrincipleComponents <- as.data.table(melt(PrincipleComponentsTrimmed,
                                                id.var = c("RalId","In2LtStatus"),
                                                variable.name = "PC", value.name="ComponentValues"))

AttemptV1 <- lm(V1 ~ In2LtStatus, data = PrincipleComponentsTrimmed)

A <- foreach(i=1:230) %do% {
   #i = 1

   #X <- lm(ComponentValues ~ In2LtStatus, data = MeltedPrincipleComponents[PC == paste("V", i, sep = "")]) #intercept = INV data
   #out = data.frame(
   #    PC = paste("PC", i, sep = ""),
   #    INV = X$coefficients[1],
   #    ST = X$coefficients[3]
   #)

   X <- lm(ComponentValues ~ In2LtStatus - 1, data = MeltedPrincipleComponents[PC == paste("V", i, sep = "")]) ### no intercept model
   out <- data.table(PC=i,
                     est=summary(X)$coef[,1],
                     se=summary(X)$coef[,2],
                     t=summary(X)$coef[,3],
                     p=summary(X)$coef[,4],
                     genotype=c("INV", "HET", "STD")) ### <- This way out structuring the output will make your life easier down the road. The only thing to double check is that the output of the summary(X) table is in the same order as you list here.

   return(out)
}

Aframe = rbindlist(A)
MeltedAFrame <- as.data.table(melt(Aframe, id.variables = c("PC"), variable.name = "InversionStatus", value.name = "coefficients"))
PCAContributionsGraph <- ggplot(data=MeltedAFrame, aes(x=PC, y=coefficients, group=InversionStatus, color = InversionStatus)) +
   geom_line() +
   geom_point()
print(PCAContributionsGraph + ggtitle("PC Coefficients from INV & ST Lines"))

'~~~~~ Graphing Singificance of Hed and STD terms ~~~~~~~~~~~~~~~~~~~'
#graph finding the significance of the hed and std terms
B <- foreach(i=1:100,.combine = "rbind") %do% {
   #i = 1
   X <- lm(ComponentValues ~ In2LtStatus, data = MeltedPrincipleComponents[PC == paste("V", i, sep = "")]) #intercept = INV data
   Y <- summary(X)$coefficients
   out = data.frame(
      t = Y[,3],
      p = Y[,4],
      term = c("Int", "het", "std"),
      id = i
   )
   return(out)
}
B <- as.data.table(B)

PCASTandINVGraph <- ggplot(B[term != "Int"], aes(x=id, y=-log10(p), group=term, color = term)) +
   geom_line() +
   geom_point() +
   geom_hline(yintercept = 0.89) + #significance of 0.05060999
   geom_hline(yintercept = 1.12) #significance of -0.04921802

print(PCASTandINVGraph + ggtitle("Significance of 'Het' and 'Std' Terms"))



'~~~ ST Bootstrap comparison with INV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
library(foreach)

#eliminated ralID from INVData:
INVDataNoRal <- subset(INVData, select = -c(ral_id))

#Shearing INVCor:

INVCorEx <- INVCor #mean = 0.05146364
INVCorEx[upper.tri(INVCorEx)] <- NA #mean = NA
diag(INVCorEx) <- NA #mean = NA
INVCorSheared <- na.omit(expand.grid(INVCorEx))

#crating a clipped ST dataset:
INVPheno <- as.vector(colnames(INVDataNoRal)) #77 Phenotypes
STDataMatchedINV <- STData[, ..INVPheno] #eliminates all phenotypes other than those mentioned in INVData

#Using foreach to get hist of ST Bootstrapped recombinations
STBootstrapINVPheno <- foreach(i=1:100) %do% {
   set.seed(i)
   nLines <- dim(STDataMatchedINV)[1]
   Experiment <- cor(STDataMatchedINV[sample(nLines, 15, replace=F),])
   Experiment[upper.tri(Experiment)] <- NA
   diag(Experiment) <- NA #removes diagonal 1 cor
   out <- data.table(
      expand.grid(Experiment)[,1],
      bootstrap.id = i,  #all outputs are labeled at bootstrap.id = i #use STBootstrapINVPheno as properly counted output
      p1 = rep(dimnames(Experiment)[[1]], (dim(Experiment)[1])),
      p2 = rep(dimnames(Experiment)[[1]], each=(dim(Experiment)[1]))
   )
  return( na.omit(out))
} #'why are there NAs?'

#creating an rbindlist from loop output
o= rbindlist(STBootstrapINVPheno) #'ral id' column included
means = o #tiddyverse to find means based on bootstrap groups group(i) command

#Boxplot Comparing Mean Averages Between Bootstraps and INV Data
means = o %>%
   group_by(bootstrap.id) %>%
   summarise(Mean = mean(V1))
mean.observed = mean(INVCorSheared$Var1)

ggplot()  +
   geom_boxplot(data = (means), aes(bootstrap.id, Mean)) +
   geom_point(aes(x = 50, y = mean.observed, size=5, col = rgb(0,0,1,0.4))) +
   ggtitle("INV/INV and ST/ST Bootstraps Correlation")
       #Warning message:
          #Continuous x aesthetic -- did you forget aes(group=...)?

   #means = as.numeric(means)
   #means = o %>%
      #group_by(bootstrap.id) %>%
      #summarise(mean.v1 = mean(o$V1)) %>%
      #mutate(bootstrap.id = "bootstrap.id") #dunno what this line did
   #mean.observed = mean(INVCorSheared$Var1)

   #ggplot()  +
     #geom_boxplot(data = (means), aes(bootstrap.id, mean.v1)) + #change the perm.id part of this
     #geom_point(aes(x = 1, y = mean.observed, size=5, col = rgb(0,0,1,0.4)))

#GGPlot Correlation Densities Bootstraps + INV Data
ggplot() +
   geom_density(data = o, aes(x = V1, group=as.factor(bootstrap.id), alpha = 0.7)) +
   geom_density(data = INVCorSheared, aes(x = Var1, col = rgb(0,0,1,0.4))) +
   ggtitle("INV/INV and ST/ST Bootstrap Correlation Densities")
   #, size = 20

#overlapping the histograms
hist(INVCorSheared, main = 'distribution of INV values',
     breaks = seq(-1.0, 1.0, 0.05)) #why isn't this breaking the histogram, or adding the title, as commanded?
histST <- lapply(STPermutationsINVPheno, hist, plot = FALSE)
histINV <- hist(INVCorSheared, plot = FALSE)
plot (histST, col = rgb(1,0,0,0.4), xlab = 'Correlation',freq = FALSE, main = 'Comparative Histogram')


histINV <- hist(INVCorSheared, breaks=10)

'~~~ Comparing one CorMatrix to another ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

STCorStripped <- fread("/Users/dgg/Desktop/STLinesComparisonMatrix.csv")
INVCorStripped <-fread("/Users/dgg/Desktop/INVLinesComparisonMatrix.csv")

cor(c(as.matrix(STCorStripped)), c(as.matrix(INVCorStripped))) #sums both matrices into vectors, compares those vectors
   0.4537692

#using Conanical Correlation Analysis (CCA) <- aims to reduce dimensionality similarly to PCA, but for two high dimensional datasets
#tries to find directions or projections that account for most of co-variance between two datasets
library(CCA)

cc_results <- cancor(STCorStripped, INVCorStripped)
cc_results$cor


'Not sure - scrapwork'

attemps <- as.data.frame(STPermutationsINVPheno)
write_xlsx(attemps,"/Users/dgg/Desktop/experiment1ST.xlsx")

INVCor
barchart(INVCor)
experiment1 <- as.data.frame(INVCor)
write_xlsx(experiment1,"/Users/dgg/Desktop/experiment1.xlsx")
