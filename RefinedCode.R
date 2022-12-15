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

### load
setwd("~/Desktop")
pheno <- readRDS("/Users/dgg/Desktop/wideform.phenotypedata (1).RDS")
write_xlsx(pheno,"/Users/dgg/Desktop\\AllPhenoNoMeans.xlsx")

'~~~Take data, subtract columns, calculate column means~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  
STwide <-fread("/Users/dgg/Desktop/ST-Lines.csv") #231 columns, 158 rows; threshold for column inclusion is 79 (50%) non-NA values in rows
STwide <- subset(STwide, select = -c(V232))
x <- as.vector(colSums(is.na(STwide)) < 79)
STnarrow <- STwide[, ..x]
STData <- STnarrow %>% 
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) #maintained 221 columns

Mixedwide <- fread("/Users/dgg/Desktop/INVST-Lines.csv") #231 columns, 25 rows; threshold for column inclusion is 13 (~50%) non-NA values in rows
Mixedwide <- subset(Mixedwide, select = -c(V232))
y <- as.vector(colSums(is.na(Mixedwide)) < 13)
STnarrow <- STwide[, ..y]
MixedData <- STnarrow %>% 
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) #maintained 213 columns


INVwide <- fread("/Users/dgg/Desktop/INV-Lines.csv") #231 columns, 18 rows; threshold for column inclusion is 9 (50%) non-NA values in rows
INVwide <- subset(INVwide, select = -c(V232))
z <- as.vector(colSums(is.na(INVwide)) < 79)
INVnarrow <- INVwide[, ..z]
INVData <- INVnarrow %>% 
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) #maintained 231 columns

'~~~Correlations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
STCor <- cor(STData)
MixedCor <- cor(MixedData)
INVCor <- cor(INVData) #In cor(INVData) : the standard deviation is zero

'~~~Heatmaps with correlation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
melted_STCor <- melt(STCor)
ggplot(data = melted_STCor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

melted_MixedCor <- melt(MixedCor)
ggplot(data = melted_MixedCor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

melted_INVCor <- melt(INVCor)
ggplot(data = melted_INVCor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

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

STDataPCA = PCA(STData, scale.unit=TRUE, ncp=5, graph=T)
res.ST <- PCA(STData, graph = FALSE)
fviz_contrib(res.ST, choice = "var", axes = 1, top = 10) #~ explains the variation from each phenotype

  #trimming the DDT studies 

MixedDataPCA = PCA(MixedData, scale.unit=TRUE, ncp=5, graph=T)
res.Mixed <- PCA(MixedData, graph = FALSE)
fviz_contrib(res.Mixed, choice = "var", axes = 1, top = 10) #~ explains the variation from each phenotype

INVDataPCA = PCA(INVData, scale.unit=TRUE, ncp=5, graph=T)
res.INV <- PCA(INVData, graph = FALSE)
fviz_contrib(res.INV, choice = "var", axes = 1, top = 10) #~ explains the variation from each phenotype

'~~~Creating Histogram~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
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

'~~~Non-linear dependencies~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' <- Is there  way to do this for a full dataset?
library(devtools) 
library(ProcessMiner/nlcor)
library(nlcor)
c <- nlcor(STData$`1Week-AcclimationSurvival_NegativeSixDegreesC_F`, STData$ActivityLevel_InducedActivity_F, plt = T)

'~~~Mutual Independence (MI) index~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
install.packages('infotheo')
library(infotheo)
MISTData <- mutinformation(STData, method="emp") #mutual information measures general dependence (including non-linear relations)








