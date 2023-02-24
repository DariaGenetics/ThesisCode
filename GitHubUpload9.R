#Note - project in 'ThesisCode' repository
setwd("/Users/dgg/Desktop")
### libraries
install.packages(...)
library(data.table)
library(ggplot2)
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
library(shinyaframe)
library(ggplot2)
library(foreach)

#PCA analysis
library(effectsize)
library(boot)
library(doMC)
registerDoMC(4)
library(patchwork)
library(missMDA)
library(factoextra)

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~In2Lt-Status~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Seasonal-Clines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

### load
setwd("~/Desktop")
pheno <- readRDS("/Users/dgg/Desktop/wideform.phenotypedata (1).RDS")
write_xlsx(pheno,"/Users/dgg/Desktop\\AllPhenoNoMeans.xlsx")

'~~~Check for Duplicates~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#find duplicated columns
duplicated_cols <- duplicated(as.data.frame(pheno))

'~~~Take data, subtract columns, calculate column means~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
Phenowide <- readRDS("/Users/dgg/Desktop/wideform.phenotypedata (1).RDS") #201 rows; threshold for column inclusion is 100 (50%) non-NA values in rows
c <- subset(Phenowide, select = -c(ral_id))
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
library(modMax)
library(igraph)

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

'~~~~~~Intercept and Inversion Status Graphs~~~~~~~~~~~~~~~~~~~~~~~~'
#incercept and inversions status flux graphs( i = 1:4, 1:10, 1:20, 1:40)
write_xlsx(PrincipleComponents,"/Users/dgg/Desktop/PrincipleComponents.xlsx") 

  'for all pheno'
PrincipleComponentsTrimmed <- fread("/Users/dgg/Desktop/PrincipleComponentsWOLastFourRows.csv") #extra 3 rows removed 
  AttemptV1 <- lm(V1 ~ In2LtStatus, data = PrincipleComponentsTrimmed)
MeltedPrincipleComponents <- as.data.table(melt(PrincipleComponentsTrimmed, id.variables = c("RalId","In2LtStatus"), variable.name = "PC", value.name="ComponentValues")) #creating really weird melting

  'for 78 phenotypes with 33% NA threshold'
PrincipleComponentsTrimmed78 <- fread("/Users/dgg/Desktop/PrincipleComponents78WOLastFourRows.csv") #extra 3 rows removed 
  Attempt78V1 <- lm(V1 ~ In2LtStatus, data = PrincipleComponentsTrimmed78)
MeltedPrincipleComponents78 <- as.data.table(melt(PrincipleComponentsTrimmed78, id.variables = c("RalId","In2LtStatus"), variable.name = "PC", value.name="ComponentValues")) #creating really weird melting
  
  
A <- foreach(i=1:20) %do% {
   #i = 1
   X <- lm(ComponentValues ~ In2LtStatus, data = MeltedPrincipleComponents[PC == paste("V", i, sep = "")]) #intercept = INV data 
   out = data.frame(
      PC = paste("PC", i, sep = ""),
      INV = X$coefficients[1],
      ST = X$coefficients[3]
   )
   return(out)
}


Aframe = rbindlist(A)
MeltedAFrame <- as.data.table(melt(Aframe, id.variables = c("PC"), variable.name = "InversionStatus", value.name = "coefficients"))
PCAContributionsGraph <- ggplot(data=MeltedAFrame, aes(x=PC, y=coefficients, group=InversionStatus, color = InversionStatus)) +
   geom_line() +
   geom_point()
print(PCAContributionsGraph + ggtitle("PC Coefficients from INV & ST Lines"))


#Non-intercept model of inversion effects on PC 
  'for all pheno'
AA <- foreach(i=1:100) %do% {
   
   XX <- lm(ComponentValues ~ In2LtStatus - 1, data = MeltedPrincipleComponents[PC == paste("V", i, sep = "")]) ### no intercept model
   out <- data.table(PC=i,
                     est=summary(XX)$coef[,1],
                     se=summary(XX)$coef[,2],
                     t=summary(XX)$coef[,3],
                     p=summary(XX)$coef[,4],
                     genotype=c("INV", "HET", "STD")) ### <- This way out structuring the output will make your life easier down the road. The only thing to double check is that the output of the summary(XX) table is in the same order as you list here.
   
   return(out)
}

AAframe = rbindlist(AA)
AAframeClipped <- AAframe[,-2:-4] #trimmed so that the extra data 'es,' 'est' aren't built into the graph
AAframeClipped$PC <- as.character(AAframeClipped$PC) 
 

MeltedAAClippedFrame <- as.data.table(melt(AAframeClipped, id.variables = c("PC"), variable.name = "genotype", value.name = "p")) #inversion status, coefficients off
MeltedAAClippedFrame2 <- MeltedAAClippedFrame[,-3] #gets rid of repeat 'p' column

PCAInversionsGraph <- ggplot(data=MeltedAAClippedFrame2, aes(x=PC, y=p, group=genotype, color = genotype)) + 
   geom_line() +
   geom_point() +
   geom_line(y = 0.05)
print(PCAInversionsGraph + ggtitle("PC P-Values by Inversion Status"))


'~~~~~ using Emmeans & Magrittr Packages ~~~~~~~~~~~~~~~~~~~'
install.packages("https://CRAN.R-project.org/package=emmeans") #package ‘https://CRAN.R-project.org/package=emmeans’ is not available (for R version 4.0.2)
  R.Version() #$major 4, $minor 0.2 - downloaded 4.2.2? (https://cran.r-project.org/src/base/R-4/)

library(emmeans) 
library(magrittr) #complementary package to Emmeans, 
   #response variable: p value 
   #factor of interest: inversion status 
   #treatment (INV) vs. control (ST): trt.vs.ctrl
  
  'possible pathway 1:'
    AllData.lm <- emmeans(AAFrame, "genotype") # "source"? # https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html#pairwise
    #AllData.lm <- emmeans(AAframe, spec = "genotype") 
    pairs(AllData.lm)

  'possible pathway 2:'
    AllData.lm <- emmeans(AAFrame, specs = trt.vs.ctrlk ~ genotype) #treatment vs. control (INV vs. ST) comparison, https://aosmith.rbind.io/2019/04/15/custom-contrasts-emmeans/#reasons-for-custom-comparisons



'~~~~~ Graphing Singificance of Hed and STD terms ~~~~~~~~~~~~~~~~~~~'
#graph finding the significance of the hed and std terms
  'for all pheno'
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

PCASTandINVGraph <- ggplot(B [term != "Int"], aes(x=id, y=-log10(p), group = term, color = term)) +
   geom_line() + 
   geom_point() + 
   geom_hline(yintercept = 0.89) + #significance of 0.05060999
   geom_hline(yintercept = 1.12) #significance of -0.04921802

print(PCASTandINVGraph + ggtitle("Significance of 'Het' and 'Std' Terms"))

  'for 78 phenotypes with 33% NA threshold'
B78 <- foreach(i=1:79,.combine = "rbind") %do% {
  #i = 1
  X <- lm(ComponentValues ~ In2LtStatus, data = MeltedPrincipleComponents78[PC == paste("V", i, sep = "")]) #intercept = INV data 
  Y <- summary(X)$coefficients
  out = data.frame( 
    t = Y[,3],
    p = Y[,4],
    term = c("Int", "het", "std"),
    id = i
  )
  return(out)
}
B78 <- as.data.table(B78)
PCASTandINVGraph78 <- ggplot(B78 [term != "Int"], aes(x=id, y=-log10(p), group=term, color = term)) +
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = 0.89) + #significance of 0.05060999
  geom_hline(yintercept = 1.12) #significance of -0.04921802
print(PCASTandINVGraph78 + ggtitle("Significance of 'Het' and 'Std' Terms - 78 Pheno"))


'~~~ ST Bootstrap comparison with INV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

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


'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~All Inversion Load & Clean~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'making a principle component graphs'
### load data
setwd("/Users/dgg/Desktop")
AllPhenoWAllInversions <- fread("AllPheno-AllInversions.csv", fill = TRUE) 

### clean up
RalID <- as.vector(unlist(AllPhenoWAllInversions$ral_id))
  
  'finding inversion counts'
  #seasonal inversion 
  In2LT_Status <- as.vector(unlist(AllPhenoWAllInversions$In2Lt)) #INV: 18, INV/ST: 25, ST: 158
    table(In2LT_Status) #<- finds INV, INV/ST, ST count for each 
  In2RNS_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RNS)) #INV: 6, INV/ST: 10, ST: 185
  In3RP_Status <- as.vector(unlist(AllPhenoWAllInversions$In3RP)) #INV: 4, INV/ST: 6, ST: 191
  
  #laditudinal cline inversions 
  In3LP_Status <- as.vector(unlist(AllPhenoWAllInversions$In3LP)) #INV: 2, INV/ST: 2, ST: 197
  In3RMo_Status <- as.vector(unlist(AllPhenoWAllInversions$In3RMo)) #INV: 17, INV/ST: 9, ST: 175
  'not sure which is In(3R)Payne - is that also In3RP, in which case expressing both seasonal and laditudinal clines?'
  
  #misc other inversions
  In2RY1_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY1))  #INV: 0, INV/ST: 1, ST: 200
  In2RY2_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY2))  #INV: 0, INV/ST: 1, ST: 200
  In2RY3_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY3))  #INV: 0, INV/ST: 1, ST: 200
  In2RY4_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY4))  #INV: 0, INV/ST: 1, ST: 200
  In2RY5_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY5))  #INV: 0, INV/ST: 1, ST: 200
  In2RY6_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY6))  #INV: 0, INV/ST: 0, ST: 201
  In2RY7_Status <- as.vector(unlist(AllPhenoWAllInversions$In2RY7))  #INV: 0, INV/ST: 0, ST: 201
  
  In3LM_Status <- as.vector(unlist(AllPhenoWAllInversions$In3LM)) #INV: 0, INV/ST: 1, ST: 200
  In3LY_Status <- as.vector(unlist(AllPhenoWAllInversions$In3LY)) #INV: 0, INV/ST: 1, ST: 200
  In3RC_Status <- as.vector(unlist(AllPhenoWAllInversions$In3RC)) #INV: 0, INV/ST: 2, ST: 199
  In3RK_Status <- as.vector(unlist(AllPhenoWAllInversions$In3RK)) #INV: 3, INV/ST: 10, ST: 188 'overdominance'
  
  UniqueInversions <- AllPhenoWAllInversions[,-c(1:231)]
  UniqueInversionStatusCountTable <- as.data.table(sapply(UniqueInversions, table))
  UniqueInversionStatusCountTable <- UniqueInversionStatusCountTable %>% rename("Status-In2RY1"="In2RY1.V1",
                                                                                "Status-In2RY2"="In2RY2.V1",
                                                                                "Status-In2RY3"="In2RY3.V1",
                                                                                "Status-In2RY4"="In2RY4.V1",
                                                                                "Status-In2RY5"="In2RY5.V1",
                                                                                "Status-In2RY6"="In2RY6.V1",
                                                                                "Status-In2RY7"="In2RY7.V1",
                                                                                "Status-In3LP"="In3LP.V1",
                                                                                "Status-In3LM"="In3LM.V1",
                                                                                "Status-In3LY"="In3LY.V1",
                                                                                "Status-In3RMo"="In3RMo.V1",
                                                                                "Status-In3RC"="In3RC.V1",
                                                                                "Status-In2RNS"="In2RNS.V1",
                                                                                "Status-In3RP"="In3RP.V1",
                                                                                "Status-In3RK"="In3RK.V1",
                                                                                "Status-In2Lt"="In2Lt.V1", #end of status renaming
                                                                                "Count-In2RY1"="In2RY1.N",
                                                                                "Count-In2RY2"="In2RY2.N",
                                                                                "Count-In2RY3"="In2RY3.N",
                                                                                "Count-In2RY4"="In2RY4.N",
                                                                                "Count-In2RY5"="In2RY5.N",
                                                                                "Count-In2RY6"="In2RY6.N",
                                                                                "Count-In2RY7"="In2RY7.N",
                                                                                "Count-In3LP"="In3LP.N",
                                                                                "Count-In3LM"="In3LM.N",
                                                                                "Count-In3LY"="In3LY.N",
                                                                                "Count-In3RMo"="In3RMo.N",
                                                                                "Count-In3RC"="In3RC.N",
                                                                                "Count-In2RNS"="In2RNS.N",
                                                                                "Count-In3RP"="In3RP.N",
                                                                                "Count-In3RK"="In3RK.N",
                                                                                "Count-In2Lt"="In2Lt.N") #end of count renaming
  #write_xlsx(UniqueInversionStatusCountTable,"/Users/dgg/Desktop\\UniqueInversionStatusCountTable.xlsx")
  
  'PCA set-up'
  AllPhenoWOAllInversions <- AllPhenoWAllInversions[, -c(232:247)]
  AllPhenoWOAllInversions <- subset(AllPhenoWOAllInversions, select = -c(ral_id))
  AllPhenoWOAllInversions <- as.matrix(AllPhenoWOAllInversions)
  dimnames(AllPhenoWOAllInversions)[[1]] <- paste("RAL_", RalID, sep="")
  
  #### impute
  AllPheno_PCA = imputePCA(AllPhenoWOAllInversions, scale.unit=TRUE, ncp=5, graph=T)
  
  ###run PCA
  AllPheno_PCA_PCs <- PCA(AllPheno_PCA, ncp=250)
  dim(AllPheno_PCA_PCs$ind$coord)
  str(AllPheno_PCA_PCs$ind$coord)
  
  ### collect output and merge with inversion status
  pc_outAllPheno <- as.data.table(AllPheno_PCA_PCs$ind$coord)
  pc_outAllPheno[,ral_id:=row.names(AllPheno_PCA_PCs$ind$coord)]
  
  
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Seasonal-Clines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~In2Lt, In2RNS, In3RP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  'In2Lt~~~~~~~~~~~~~~'
  ### merge data
  inv_In2Lt <- AllPhenoWAllInversions[,c("ral_id", "In2Lt"), with=F]
  inv_In2Lt[,ral_id:=paste("RAL_", ral_id, sep="")]
  pc_In2Lt <- merge(pc_outAllPheno, inv_In2Lt, by="ral_id")
  
  MeltedPCs_In2Lt <- as.data.table(melt(pc_In2Lt,
                                        id.var = c("ral_id","In2Lt"),
                                        variable.name = "PC", value.name="ComponentValues"))
  
  
  
  ### iterate through pcs
  A_In2Lt <- foreach(i=unique(MeltedPCs_In2Lt$PC), .combine="rbind") %dopar% {
    message(i)
    #i = unique(MeltedPrincipleComponents$PC)[100]
    
    Xa_In2Lt <- lm(ComponentValues ~ In2Lt, data = MeltedPCs_In2Lt[PC == i])
    
    foo_In2Lt <- boot(MeltedPCs_In2Lt[PC == i],
                function(data, indices) summary(lm(ComponentValues~In2Lt, data[indices,]))$r.squared,
                R=1000)
    
    
    out_In2Lt <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                      lci=quantile(foo_In2Lt$t, c(0.025)),
                      uci=quantile(foo_In2Lt$t, c(0.975)),
                      r2=summary(Xa_In2Lt)$r.squared,
                      p=summary(aov(Xa_In2Lt))[[1]][1,5],
                      pve=as.numeric(as.data.table(AllPheno_PCA_PCs$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
    
    return(out_In2Lt)
  }
  A_In2Lt[,pa:=p.adjust(p)]
  
  r2_plot_In2Lt <-
    ggplot(data=A_In2Lt, aes(x=PC, y=r2)) +
    geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
    geom_point() +
    geom_point(data=A_In2Lt[p<.05], aes(x=PC, y=r2), color="red") +
    geom_point(data=A_In2Lt[pa<.1], aes(x=PC, y=r2), color="green") +
    ggtitle("Principle Components - In2LT") + 
    theme(plot.title = element_text(size=12, hjust = 0.5))
  
  ### which factors contribute to Dim.8?
  loading_In2Lt <- as.data.table(AllPheno_PCA_PCs$var$coord)
  loading_In2Lt[,pheno:=row.names(AllPheno_PCA_PCs$var$coord)]
  ll_In2Lt <- melt(loading_In2Lt, id.var="pheno")
  
  ll.rank_In2Lt <- ll_In2Lt[,list(quan=rank(abs(value))/(length(value)+1), pheno=pheno, rank=rank(value)), list(variable)]
  
  setkey(ll_In2Lt, variable, pheno)
  setkey(ll.rank_In2Lt, variable, pheno)
  ll_In2Lt <- merge(ll_In2Lt, ll.rank_In2Lt)
  
  'full blue - 8'
  'full red - 4, 9, 112, 128, 137, 145, 179, 191, 197'
  
  #contributing factors
  #PC8 - Blue
  phenoLoad_In2Lt_8 <-
    ggplot(data=ll_In2Lt[variable=="Dim.8"][quan>.9]) +
    geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
    geom_text(data=ll_In2Lt[variable=="Dim.8"][quan>.9],
              aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
    xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
    ggtitle("Top 10% Influencing Phenotypes: PC 8
      p<0.05 - Bonferroni Correction") + 
    theme(plot.title = element_text(size=8, hjust = 0.5))
  
  #PC4 - Red
  phenoLoad_In2Lt_4 <-
    ggplot(data=ll_In2Lt[variable=="Dim.4"][quan>.9]) +
    geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
    geom_text(data=ll_In2Lt[variable=="Dim.4"][quan>.9],
              aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
    xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
    ggtitle("Top 10% Influencing Phenotypes: PC 4
      p<0.05 - ANOVA Test") + 
    theme(plot.title = element_text(size=8, hjust = 0.5))
  
  #PC9 - Red
  phenoLoad_In2Lt_9 <-
    ggplot(data=ll_In2Lt[variable=="Dim.9"][quan>.9]) +
    geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
    geom_text(data=ll_In2Lt[variable=="Dim.9"][quan>.9],
              aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
    xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
    ggtitle("Top 10% Influencing Phenotypes: PC 9
      p<0.05 - ANOVA Test") + 
    theme(plot.title = element_text(size=8, hjust = 0.5))
  
  ### effect plot
  pc.ag_In2Lt <- MeltedPCs_In2Lt[, list(mu=mean(ComponentValues), sd=sd(ComponentValues), n=length(ComponentValues)),
                                     list(In2Lt, PC)]
  pc.ag_In2Lt[,se:=sd/sqrt(n)]
  
  effect_In2Lt_8 <- ggplot(data=pc.ag_In2Lt[PC=="Dim.8"], aes(x=In2Lt, y=mu)) +
    geom_point() +
    geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se)) +
    ggtitle("Effect Size: PC 8
      p<0.05 - Bonferroni Correction") + 
    theme(plot.title = element_text(size=8, hjust = 0.5))
  
  effect_In2Lt_4 <- ggplot(data=pc.ag_In2Lt[PC=="Dim.4"], aes(x=In2Lt, y=mu)) +
    geom_point() +
    geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se)) +
    ggtitle("Effect Size: PC 4
      p<0.05 - ANOVA Test") + 
    theme(plot.title = element_text(size=8, hjust = 0.5))
  
  effect_In2Lt_9 <- ggplot(data=pc.ag_In2Lt[PC=="Dim.9"], aes(x=In2Lt, y=mu)) +
    geom_point() +
    geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se)) +
    ggtitle("Effect Size: PC 9
      p<0.05 - ANOVA Test") + 
    theme(plot.title = element_text(size=8, hjust = 0.5))
  
  ##create plot 
  layout <- "
  AAAA
  BBCC
  DDEE
  FFGG"
  
  r2_plot_In2Lt + phenoLoad_In2Lt_8 + effect_In2Lt_8 + phenoLoad_In2Lt_4 + effect_In2Lt_4 + phenoLoad_In2Lt_9 + effect_In2Lt_9 +
    plot_layout(design=layout) + plot_annotation(tag_levels ="A")
  
  
  'In2RNS~~~~~~~~~~~~~~'
  ### merge data
  inv_In2RNS <- AllPhenoWAllInversions[,c("ral_id", "In2RNS"), with=F]
  inv_In2RNS[,ral_id:=paste("RAL_", ral_id, sep="")]
  pc_In2RNS <- merge(pc_outAllPheno, inv_In2RNS, by="ral_id")
  
  MeltedPCs_In2RNS <- as.data.table(melt(pc_In2RNS,
                                        id.var = c("ral_id","In2RNS"),
                                        variable.name = "PC", value.name="ComponentValues"))
 
   ### iterate through pcs
  A_In2RNS <- foreach(i=unique(MeltedPCs_In2RNS$PC), .combine="rbind") %dopar% {
    message(i)
    #i = unique(MeltedPrincipleComponents$PC)[100]
    
    Xa_In2RNS <- lm(ComponentValues ~ In2RNS, data = MeltedPCs_In2RNS[PC == i])
    
    foo_In2RNS <- boot(MeltedPCs_In2RNS[PC == i],
                      function(data, indices) summary(lm(ComponentValues~In2RNS, data[indices,]))$r.squared,
                      R=1000)
    
    
    out_In2RNS <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                            lci=quantile(foo_In2RNS$t, c(0.025)),
                            uci=quantile(foo_In2RNS$t, c(0.975)),
                            r2=summary(Xa_In2RNS)$r.squared,
                            p=summary(aov(Xa_In2RNS))[[1]][1,5],
                            pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
    
    return(out_In2RNS)
  }
  A_In2RNS[,pa:=p.adjust(p)]
  
  r2_plot_In2RNS <-
    ggplot(data=A_In2RNS, aes(x=PC, y=r2)) +
    geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
    geom_point() +
    geom_point(data=A_In2RNS[p<.05], aes(x=PC, y=r2), color="red") +
    geom_point(data=A_In2RNS[pa<.1], aes(x=PC, y=r2), color="blue") +
    ggtitle("Principle Components - In2RNS") + 
    theme(plot.title = element_text(size=12, hjust = 0.5))
  
  'In3RP~~~~~~~~~~~~~~'
  ### merge data
  inv_In3RP <- AllPhenoWAllInversions[,c("ral_id", "In3RP"), with=F]
  inv_In3RP[,ral_id:=paste("RAL_", ral_id, sep="")]
  pc_In3RP <- merge(pc_outAllPheno, inv_In3RP, by="ral_id")
  
  MeltedPCs_In3RP <- as.data.table(melt(pc_In3RP,
                                         id.var = c("ral_id","In3RP"),
                                         variable.name = "PC", value.name="ComponentValues"))
  
  ### iterate through pcs
  A_In3RP <- foreach(i=unique(MeltedPCs_In3RP$PC), .combine="rbind") %dopar% {
    message(i)
    #i = unique(MeltedPrincipleComponents$PC)[100]
    
    Xa_In3RP <- lm(ComponentValues ~ In3RP, data = MeltedPCs_In3RP[PC == i])
    
    foo_In3RP <- boot(MeltedPCs_In3RP[PC == i],
                       function(data, indices) summary(lm(ComponentValues~In3RP, data[indices,]))$r.squared,
                       R=1000)
    
    
    out_In3RP <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                             lci=quantile(foo_In3RP$t, c(0.025)),
                             uci=quantile(foo_In3RP$t, c(0.975)),
                             r2=summary(Xa_In3RP)$r.squared,
                             p=summary(aov(Xa_In3RP))[[1]][1,5],
                             pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
    
    return(out_In3RP)
  }
  A_In3RP[,pa:=p.adjust(p)]
  
  r2_plot_In3RP <-
    ggplot(data=A_In3RP, aes(x=PC, y=r2)) +
    geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
    geom_point() +
    geom_point(data=A_In3RP[p<.05], aes(x=PC, y=r2), color="red") +
    geom_point(data=A_In3RP[pa<.1], aes(x=PC, y=r2), color="blue") +
    ggtitle("Principle Components - In3RP") + 
    theme(plot.title = element_text(size=12, hjust = 0.5))
  
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~Latitudinal-Clines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~In3LP, In3RMo, In3RPayne(?)~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

'In3LP~~~~~~~~~~~~~~'
### merge data
inv_In3LP <- AllPhenoWAllInversions[,c("ral_id", "In3LP"), with=F]
inv_In3LP[,ral_id:=paste("RAL_", ral_id, sep="")]
pc_In3LP <- merge(pc_outAllPheno, inv_In3LP, by="ral_id")

MeltedPCs_In3LP <- as.data.table(melt(pc_In3LP,
                                      id.var = c("ral_id","In3LP"),
                                      variable.name = "PC", value.name="ComponentValues"))

### iterate through pcs
A_In3LP <- foreach(i=unique(MeltedPCs_In3LP$PC), .combine="rbind") %dopar% {
  message(i)
  #i = unique(MeltedPrincipleComponents$PC)[100]
  
  Xa_In3LP <- lm(ComponentValues ~ In3LP, data = MeltedPCs_In3LP[PC == i])
  
  foo_In3LP <- boot(MeltedPCs_In3LP[PC == i],
                    function(data, indices) summary(lm(ComponentValues~In3LP, data[indices,]))$r.squared,
                    R=1000)
  
  
  out_In3LP <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                          lci=quantile(foo_In3LP$t, c(0.025)),
                          uci=quantile(foo_In3LP$t, c(0.975)),
                          r2=summary(Xa_In3LP)$r.squared,
                          p=summary(aov(Xa_In3LP))[[1]][1,5],
                          pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
  
  return(out_In3LP)
}
A_In3LP[,pa:=p.adjust(p)]

r2_plot_In3LP <-
  ggplot(data=A_In3LP, aes(x=PC, y=r2)) +
  geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
  geom_point() +
  geom_point(data=A_In3LP[p<.05], aes(x=PC, y=r2), color="red") +
  geom_point(data=A_In3LP[pa<.1], aes(x=PC, y=r2), color="blue") +
  ggtitle("Principle Components - In3LP") + 
  theme(plot.title = element_text(size=12, hjust = 0.5))

'In3RMo~~~~~~~~~~~~~~'
### merge data
inv_In3RMo <- AllPhenoWAllInversions[,c("ral_id", "In3RMo"), with=F]
inv_In3RMo[,ral_id:=paste("RAL_", ral_id, sep="")]
pc_In3RMo <- merge(pc_outAllPheno, inv_In3RMo, by="ral_id")

MeltedPCs_In3RMo <- as.data.table(melt(pc_In3RMo,
                                      id.var = c("ral_id","In3RMo"),
                                      variable.name = "PC", value.name="ComponentValues"))

### iterate through pcs
A_In3RMo <- foreach(i=unique(MeltedPCs_In3RMo$PC), .combine="rbind") %dopar% {
  message(i)
  #i = unique(MeltedPrincipleComponents$PC)[100]
  
  Xa_In3RMo <- lm(ComponentValues ~ In3RMo, data = MeltedPCs_In3RMo[PC == i])
  
  foo_In3RMo <- boot(MeltedPCs_In3RMo[PC == i],
                    function(data, indices) summary(lm(ComponentValues~In3RMo, data[indices,]))$r.squared,
                    R=1000)
  
  
  out_In3RMo <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                          lci=quantile(foo_In3RMo$t, c(0.025)),
                          uci=quantile(foo_In3RMo$t, c(0.975)),
                          r2=summary(Xa_In3RMo)$r.squared,
                          p=summary(aov(Xa_In3RMo))[[1]][1,5],
                          pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
  
  return(out_In3RMo)
}
A_In3RMo[,pa:=p.adjust(p)]

r2_plot_In3RMo <-
  ggplot(data=A_In3RMo, aes(x=PC, y=r2)) +
  geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
  geom_point() +
  geom_point(data=A_In3RMo[p<.05], aes(x=PC, y=r2), color="red") +
  geom_point(data=A_In3RMo[pa<.1], aes(x=PC, y=r2), color="blue") +
  ggtitle("Principle Components - In3RMo") + 
  theme(plot.title = element_text(size=12, hjust = 0.5))

r2_plot_In3RMo

### which factors contribute to Dim.8?
loading_In3RMo <- as.data.table(AllPheno_PCA_PCs$var$coord)
loading_In3RMo[,pheno:=row.names(AllPheno_PCA_PCs$var$coord)]
ll_In3RMo <- melt(loading_In3RMo, id.var="pheno")

ll.rank_In3RMo <- ll_In3RMo[,list(quan=rank(abs(value))/(length(value)+1), pheno=pheno, rank=rank(value)), list(variable)]

setkey(ll_In3RMo, variable, pheno)
setkey(ll.rank_In3RMo, variable, pheno)
ll_In3RMo <- merge(ll_In3RMo, ll.rank_In3RMo)

#blue: dim 4
#red: 20, 22, 39, 44, 47, 171, 189, 196, 198

#contributing factors
phenoLoad_In3RMo_4 <-
  ggplot(data=ll_In3RMo[variable=="Dim.4"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll_In3RMo[variable=="Dim.4"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
  ggtitle("Top 10% Influencing Phenotypes: PC 4
      p<0.05 - Bonferroni Correction") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

phenoLoad_In3RMo_20 <-
  ggplot(data=ll_In3RMo[variable=="Dim.20"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll_In3RMo[variable=="Dim.20"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")+
  ggtitle("Top 10% Influencing Phenotypes: PC 20
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

phenoLoad_In3RMo_22 <-
  ggplot(data=ll_In3RMo[variable=="Dim.22"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll_In3RMo[variable=="Dim.22"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")+
  ggtitle("Top 10% Influencing Phenotypes: PC 22
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

#effect plots
pc.ag_In3RMo <- MeltedPCs_In3RMo[, list(mu=mean(ComponentValues), sd=sd(ComponentValues), n=length(ComponentValues)),
                                   list(In3RMo, PC)]
pc.ag_In3RMo[,se:=sd/sqrt(n)]

effect_plot_In3RMo_4 <- ggplot(data=pc.ag_In3RMo[PC=="Dim.4"], aes(x=In3RMo, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In3RMo, xend=In3RMo, y=mu-2*se, yend=mu+2*se)) +
  ggtitle("Effect Size: PC 4
      p<0.05 - Bonferroni Correction") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

effect_plot_In3RMo_20 <- ggplot(data=pc.ag_In3RMo[PC=="Dim.20"], aes(x=In3RMo, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In3RMo, xend=In3RMo, y=mu-2*se, yend=mu+2*se)) +
  ggtitle("Effect Size: PC 20
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

effect_plot_In3RMo_22 <- ggplot(data=pc.ag_In3RMo[PC=="Dim.22"], aes(x=In3RMo, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In3RMo, xend=In3RMo, y=mu-2*se, yend=mu+2*se)) +
  ggtitle("Effect Size: PC 22
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

layout <- "
  AAAA
  BBCC
  DDEE
  FFGG"

r2_plot_In3RMo + phenoLoad_In3RMo_4 + effect_plot_In3RMo_4 + 
  phenoLoad_In3RMo_20 + effect_plot_In3RMo_20 + 
  phenoLoad_In3RMo_22 + effect_plot_In3RMo_22 +
  plot_layout(design=layout) + plot_annotation(tag_levels ="A")

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~Miscallenous Inversions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~In2RY1-7, In3LM, In3LY, In3RC, In3RK (Overdominance)~~~~~~~~~~~~~~~~~~~~~~'

'In2RY1~~~~~~~~~~~~~'
### merge data
inv_In2RY1 <- AllPhenoWAllInversions[,c("ral_id", "In2RY1"), with=F]
inv_In2RY1[,ral_id:=paste("RAL_", ral_id, sep="")]
pc_In2RY1 <- merge(pc_outAllPheno, inv_In2RY1, by="ral_id")

MeltedPCs_In2RY1 <- as.data.table(melt(pc_In2RY1,
                                      id.var = c("ral_id","In2RY1"),
                                      variable.name = "PC", value.name="ComponentValues"))

### iterate through pcs
A_In2RY1 <- foreach(i=unique(MeltedPCs_In2RY1$PC), .combine="rbind") %dopar% {
  message(i)
  #i = unique(MeltedPrincipleComponents$PC)[100]
  
  Xa_In2RY1 <- lm(ComponentValues ~ In2RY1, data = MeltedPCs_In2RY1[PC == i])
  
  foo_In2RY1 <- boot(MeltedPCs_In2RY1[PC == i],
                    function(data, indices) summary(lm(ComponentValues~In2RY1, data[indices,]))$r.squared,
                    R=1000)
  
  
  out_In2RY1 <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                          lci=quantile(foo_In2RY1$t, c(0.025)),
                          uci=quantile(foo_In2RY1$t, c(0.975)),
                          r2=summary(Xa_In2RY1)$r.squared,
                          p=summary(aov(Xa_In2RY1))[[1]][1,5],
                          pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
  
  return(out_In2RY1)
}
  #Error in { : 
    #task 1 failed - "contrasts can be applied only to factors with 2 or more levels"

'In3RK~~~~~~~~~~~~~~'
### merge data
inv_In3RK <- AllPhenoWAllInversions[,c("ral_id", "In3RK"), with=F]
inv_In3RK[,ral_id:=paste("RAL_", ral_id, sep="")]
pc_In3RK <- merge(pc_outAllPheno, inv_In3RK, by="ral_id")

MeltedPCs_In3RK <- as.data.table(melt(pc_In3RK,
                                       id.var = c("ral_id","In3RK"),
                                       variable.name = "PC", value.name="ComponentValues"))

### iterate through pcs
A_In3RK <- foreach(i=unique(MeltedPCs_In3RK$PC), .combine="rbind") %dopar% {
  message(i)
  #i = unique(MeltedPrincipleComponents$PC)[100]
  
  Xa_In3RK <- lm(ComponentValues ~ In3RK, data = MeltedPCs_In3RK[PC == i])
  
  foo_In3RK <- boot(MeltedPCs_In3RK[PC == i],
                     function(data, indices) summary(lm(ComponentValues~In3RK, data[indices,]))$r.squared,
                     R=1000)
  
  
  out_In3RK <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                           lci=quantile(foo_In3RK$t, c(0.025)),
                           uci=quantile(foo_In3RK$t, c(0.975)),
                           r2=summary(Xa_In3RK)$r.squared,
                           p=summary(aov(Xa_In3RK))[[1]][1,5],
                           pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
  
  return(out_In3RK)
}
A_In3RK[,pa:=p.adjust(p)]

r2_plot_In3RK <-
  ggplot(data=A_In3RK, aes(x=PC, y=r2)) +
  geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
  geom_point() +
  geom_point(data=A_In3RK[p<.05], aes(x=PC, y=r2), color="red") +
  geom_point(data=A_In3RK[pa<.1], aes(x=PC, y=r2), color="blue") +
  ggtitle("Principle Components - In3RK") + 
  theme(plot.title = element_text(size=12, hjust = 0.5))

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~Comparing Bonferroni-Corrected Significant Principle Components~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

PhenoLoadNames_In2Lt_8=ll_In2Lt[variable=="Dim.8"][quan>.9]
PhenoLoadNames_In2Lt_8 <- PhenoLoadNames_In2Lt_8$pheno

phenoLoadNames_In3RMo_4 = ll_In3RMo[variable=="Dim.4"][quan>.9]
phenoLoadNames_In3RMo_4 <- phenoLoadNames_In3RMo_4$pheno

Bonferroni_Significant_Pheno <- data.frame(PhenoLoadNames_In2Lt_8, phenoLoadNames_In3RMo_4)
Bonferroni_Significant_Pheno_Graph <- barplot(prop.table(table(Bonferroni_Significant_Pheno)), xlab = "phenotype", ylab = "frequency", main = "Phenotype Frequency in Bonferroni-Corrected Significant Principle Components")

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~Comparing Bonferroni/ANOVA Significant Principle Components~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  #ammasing influencing phenotypes 
  
  PhenoLoadNames_In2Lt_8=ll_In2Lt[variable=="Dim.8"][quan>.9]
  PhenoLoadNames_In2Lt_8 <- PhenoLoadNames_In2Lt_8$pheno
  PhenoLoadNames_In2Lt_8DT <- as.data.frame(PhenoLoadNames_In2Lt_8)
  
  
  PhenoLoadNames_In3RMo_4=ll_In3RMo[variable=="Dim.4"][quan>.9]
  PhenoLoadNames_In3RMo_4 <- PhenoLoadNames_In3RMo_4$pheno
  PhenoLoadNames_In3RMo_4 <- as.data.frame(PhenoLoadNames_In3RMo_4)
  
  CombinedPheno <- as.data.frame(c(PhenoLoadNames_In2Lt_8DT, PhenoLoadNames_In3RMo_4))
  
  write_xlsx(CombinedPheno,"/Users/dgg/Desktop\\PhenoForGWAS.xlsx")
  
  
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~GWAS Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
#load in data - INVERSION 3RMo PHENO
TwoHeptanoneAttraction_standard_M <- fread("2-HeptanoneAttraction_standard_male.nogrms.original.txt",
                                                 stringsAsFactor = FALSE)
TwoHeptanoneAttraction_standard_F <- fread("2-HeptanoneAttraction_standard_female.nogrms.original.txt",
                                                 stringsAsFactor = FALSE)
OneHexanolAttraction_standard_F <- fread("1-HexanolAttraction_standard_female.nogrms.original.txt",
                                               stringsAsFactor = FALSE) 
CentroidSize_Standard_F <- fread("CentroidSize_standard_F.nogrms.jan.txt",
                                      stringsAsFactor = FALSE)
CentroidSize_Standard_M <- fread("CentroidSize_standard_M.nogrms.jan.txt",
                                      stringsAsFactor = FALSE)
HexanalAttraction_standard_M <- fread("HexanalAttraction_standard_male.nogrms.original.txt",
                                           stringsAsFactor = FALSE)
HexanalAttraction_standard_F <- fread("HexanalAttraction_standard_female.nogrms.original.txt",
                                           stringsAsFactor = FALSE)
IntensityOfResponse_standard_F <- fread("IntensityOfResponse_standard_F.nogrms.original.txt",
                                             stringsAsFactor = FALSE)
InternocularDistance_standard_F <- fread("InternocularDistance_standard_F.nogrms.jan.txt",
                                              stringsAsFactor = FALSE)
InternocularDistance_standard_M <- fread("InternocularDistance_standard_M.nogrms.jan.txt",
                                              stringsAsFactor = FALSE) #there's duplicates here 
LifeTimeFecundity_LowYeast_F <- fread("LifeTimeFecundity_LowYeast_female.nogrms.original.txt",
                                           stringsAsFactor = FALSE)
NegativeGeotaxis_MSB_Treatment_F <- fread("NegativeGeotaxis_MSB_Treatment_female.nogrms.original.txt",
                                               stringsAsFactor = FALSE)
PupalCaseLength_Standard_Mix <- fread("PupalCaseLength_Standard_Mix.nogrms.jan.txt",
                                           stringsAsFactor = FALSE) 
ThoraxLength_standard_F <- fread("ThoraxLength_standard_female.nogrms.original.txt",
                                      stringsAsFactor = FALSE)
TriglycerideLevels_PooledDiets_M <- fread("TriglycerideLevels_PooledDiets_male.nogrms.original.txt",
                                               stringsAsFactor = FALSE)
TriglycerideLevels_LowGlucoseDiet_M <- fread("TriglycerideLevels_LowGlucoseDiet_male.nogrms.original.txt",
                                                  stringsAsFactor = FALSE)
Weight_PooledDiets_M <- fread("Weight_PooledDiets_male.nogrms.original.txt",
                                   stringsAsFactor = FALSE)
Weight_LowGlucoseDiet_M <- fread("Weight_LowGlucoseDiet_male.nogrms.original.txt",
                                      stringsAsFactor = FALSE)
Weight_HighGlucoseDiet_M <- fread("Weight_HighGlucoseDiet_male.nogrms.original.txt",
                                       stringsAsFactor = FALSE)
WingCentroidSize_standard_M <- fread("WingCentroidSize_standard_male.nogrms.original.txt",
                                          stringsAsFactor = FALSE)
WingCentroidSize_standard_F <- fread("WingCentroidSize_standard_female.nogrms.original.txt",
                                          stringsAsFactor = FALSE)
VectorOf_3RMo_Pheno<- c("2-HeptanoneAttraction_standard_male.nogrms.original.txt", "2-HeptanoneAttraction_standard_female.nogrms.original.txt", 
                        "1-HexanolAttraction_standard_female.nogrms.original.txt", "CentroidSize_standard_F.nogrms.jan.txt", 
                        "CentroidSize_standard_M.nogrms.jan.txt", "HexanalAttraction_standard_male.nogrms.original.txt", 
                        "HexanalAttraction_standard_female.nogrms.original.txt", "IntensityOfResponse_standard_F.nogrms.original.txt", 
                        "InternocularDistance_standard_F.nogrms.jan.txt", "InternocularDistance_standard_M.nogrms.jan.txt", 
                        "LifeTimeFecundity_LowYeast_female.nogrms.original.txt", "NegativeGeotaxis_MSB_Treatment_female.nogrms.original.txt", 
                        "PupalCaseLength_Standard_Mix.nogrms.jan.txt", "ThoraxLength_standard_female.nogrms.original.txt", 
                        "TriglycerideLevels_PooledDiets_male.nogrms.original.txt", "TriglycerideLevels_LowGlucoseDiet_male.nogrms.original.txt",
                        "Weight_PooledDiets_male.nogrms.original.txt","Weight_LowGlucoseDiet_male.nogrms.original.txt",
                        "Weight_HighGlucoseDiet_male.nogrms.original.txt", "WingCentroidSize_standard_male.nogrms.original.txt",
                        "WingCentroidSize_standard_female.nogrms.original.txt")

#load in data - INVERSION 2Lt PHENO
ChillComaRecoveryTime_standard_M <- fread("ChillComaRecoveryTime_standard_male.nogrms.original.txt",
                                               stringsAsFactor = FALSE)
ChillComaRecoveryTime_standard_F <- fread("ChillComaRecoveryTime_standard_female.nogrms.original.txt",
                                               stringsAsFactor = FALSE)
ClimbingHeight_ParaquatInsecticideExposure_M <- fread("ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt",
                                                           stringsAsFactor = FALSE)
ClimbingHeight_standard_M <- fread("ClimbingHeight_standard_M.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
DaySleep_standard_F <- fread("DaySleep_standard_Female.nogrms.original.txt",
                                  stringsAsFactor = FALSE)
DayBoutNumber_standard_F <- fread("DayBoutNumber_standard_Female.nogrms.original.txt",
                                       stringsAsFactor = FALSE)
FoodIntake_standard_F <- fread("FoodIntake_standard_female.nogrms.original.txt",
                                    stringsAsFactor = FALSE)
FreeGlycerolLevels_PooledDiets_M <- fread("FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt",
                                               stringsAsFactor = FALSE)
FreeGlycerolLevels_LowGlucoseDiet_M <- fread("FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt",
                                                  stringsAsFactor = FALSE)
FreeGlycerolLevels_HighGlucoseDiet_M <- fread("FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt",
                                                   stringsAsFactor = FALSE) 
MeanElutionTime_RepeatedEthanolExposures_M <- fread("MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt",
                                                         stringsAsFactor = FALSE)
MeanElutionTime_RepeatedEthanolExposures_F <- fread("MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt",
                                                         stringsAsFactor = FALSE)
MinimumSurvivableTemperature_TwentyNineDegreesC_M <- fread("MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt",
                                                                stringsAsFactor = FALSE)
MinimumSurvivableTemperature_TwentySixDegreesC_M <- fread("MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt",
                                                               stringsAsFactor = FALSE)
MinimumSurvivableTemperature_TwentyThreeDegreesC_M <- fread("MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt",
                                                                 stringsAsFactor = FALSE)  
SolubleProteinLevels_PooledDiets_M <- fread("SolubleProteinLevels_PooledDiets_male.nogrms.original.txt",
                                                 stringsAsFactor = FALSE)
SolubleProteinLevels_LowGlucoseDiet_M <- fread("SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt",
                                                    stringsAsFactor = FALSE)  
Weight_PooledDiets_M <- fread("Weight_PooledDiets_male.nogrms.original.txt",
                                   stringsAsFactor = FALSE)
Weight_LowGlucoseDiet_M <- fread("Weight_LowGlucoseDiet_male.nogrms.original.txt",
                                      stringsAsFactor = FALSE)
Weight_HighGlucoseDiet_M <- fread("Weight_HighGlucoseDiet_male.nogrms.original.txt",
                                       stringsAsFactor = FALSE)
VectorOf_2Lt_Pheno <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "ChillComaRecoveryTime_standard_female.nogrms.original.txt",
                       "ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt", "ClimbingHeight_standard_M.nogrms.jan.txt",
                       "DaySleep_standard_Female.nogrms.original.txt", "DayBoutNumber_standard_Female.nogrms.original.txt", 
                       "FoodIntake_standard_female.nogrms.original.txt", "FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt", 
                       "FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt", "FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt", 
                       "MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt", "MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt", 
                       "MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt", "MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt", 
                       "MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt", "SolubleProteinLevels_PooledDiets_male.nogrms.original.txt", 
                       "SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt", "Weight_PooledDiets_male.nogrms.original.txt", 
                       "Weight_LowGlucoseDiet_male.nogrms.original.txt", "Weight_HighGlucoseDiet_male.nogrms.original.txt") 
  
'~~~~~~~~~~~'
setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")

library(data.table)
library(ggplot2)
library(foreach)

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#Is there enrichment in the intervsion? - EXAMPLE
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
VectorOf_2Lt_Pheno_experiment <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "ChillComaRecoveryTime_standard_female.nogrms.original.txt")

PhenoGather_EnrichmnetINV_2Lt <- foreach(i=1:length(VectorOf_2Lt_Pheno_experiment)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")

  data <- fread(VectorOf_2Lt_Pheno_experiment[i])
  
 
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_2Lt_Pheno_experiment[i]] #can use a split to take off file extensions from names
  return(data2)
}
gather_2Lt_experiment <- rbindlist(PhenoGather_EnrichmnetINV_2Lt)
#null hypothesis that inv woud have 0.001

  #for graph separating the phenotypes
ggplot(gather_2Lt_experiment, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar") +
  xlab("inv") +
  ylab("nTOP") +
  facet_wrap(~Phenotype) +
  geom_hline(yintercept=0.0010, linetype = "dashed", color ="red")

  #for graph combining the phenotypes
ggplot(gather_2Lt_experiment, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar") +
  xlab("inv") +
  ylab("nTOP") +
  geom_hline(yintercept=0.0010, linetype = "dashed", color ="red")

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#Is there enrichment in the inversion? - FOR ALL DATA EACH INVERSION
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  #INVERSION IN2LT #################################################################
'thr 0.05~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
Enrichment_FULL_INV_2Lt_thr0.05 <- foreach(i=1:length(VectorOf_2Lt_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  
  data <- fread(VectorOf_2Lt_Pheno[i])
  
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.05)), thr = c(0.05), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_2Lt_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
Enrichment_FULL_INV_2Lt_Rbound_thr0.05 <- rbindlist(Enrichment_FULL_INV_2Lt_thr0.05)

#combined phenotypes thr 0.05
ggplot_In2Lt_thr_0.05 <- ggplot(Enrichment_FULL_INV_2Lt_Rbound_thr0.05, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
  xlab("Chromosome Section") +
  ylab("nTOP") +
  geom_hline(yintercept=0.05, linetype = "dashed", color ="red") +
  facet_wrap(~thr) + 
  theme(axis.text.x=element_text(size=9)) +
  ggtitle("2Lt Enrichment, thr = 0.05")


'thr 0.01~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#null hypothesis that inv would have 0.001
Enrichment_FULL_INV_2Lt_thr0.01 <- foreach(i=1:length(VectorOf_2Lt_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  
  data <- fread(VectorOf_2Lt_Pheno[i])
  
  
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.01)), thr = c(0.01), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_2Lt_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
Enrichment_FULL_INV_2Lt_Rbound_thr0.01 <- rbindlist(Enrichment_FULL_INV_2Lt_thr0.01)

#combined phenotypes thr 0.01 
ggplot_In2Lt_thr_0.01 <- ggplot(Enrichment_FULL_INV_2Lt_Rbound_thr0.01, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
  xlab("inv") +
  ylab("nTOP") +
  geom_hline(yintercept=0.01, linetype = "dashed", color ="red") +
  facet_wrap(~thr) +
  theme(axis.text.x=element_text(size=9)) +
  ggtitle("2Lt Enrichment, thr = 0.01")


'thr 0.001~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
Enrichment_FULL_INV_2Lt_thr0.001 <- foreach(i=1:length(VectorOf_2Lt_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  
  data <- fread(VectorOf_2Lt_Pheno[i])
  
  
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_2Lt_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
Enrichment_FULL_INV_2Lt_Rbound_thr0.001 <- rbindlist(Enrichment_FULL_INV_2Lt_thr0.001)

#combined phenotypes thr 0.001 
ggplot_In2Lt_thr_0.001 <- ggplot(Enrichment_FULL_INV_2Lt_Rbound_thr0.001, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar") +
  xlab("inv") +
  ylab("nTOP") +
  geom_hline(yintercept=0.001, linetype = "dashed", color ="red") +
  facet_wrap(~thr) +
  theme(axis.text.x=element_text(size=9)) +
  ggtitle("2Lt Enrichment, thr = 0.001")


'all thr graphs, combined'
library(gridExtra)
# Arrange the plots in a grid with a 1x3 layout
grid.arrange(ggplot_In2Lt_thr_0.05, ggplot_In2Lt_thr_0.01, ggplot_In2Lt_thr_0.001, ncol = 1) 
  #+grid.text("Tag A", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = "left", gp = gpar(fontsize = 16, fontface = "bold"))

'for all pheno, separated'
  #per individual phenotype
ggplot(Enrichment_FULL_INV_2Lt_Rbound, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar") +
  xlab("inv") +
  ylab("nTOP") +
  facet_wrap(~Phenotype) +
  geom_hline(yintercept=0.0010, linetype = "dashed", color ="red") +
  ggtitle("All Phenotypes Inversion 2Lt Enrichment Compared to Null Hypothesis (nTop = 0.001)")
  



#INVERSION 3RMo #################################################################

'thr 0.05~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
Enrichment_FULL_INV_3RMo_thr0.05 <- foreach(i=1:length(VectorOf_3RMo_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  
  data <- fread(VectorOf_3RMo_Pheno[i])
  
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.05)), thr = c(0.05), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_3RMo_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
Enrichment_FULL_INV_3RMo_Rbound_thr0.05 <- rbindlist(Enrichment_FULL_INV_3RMo_thr0.05)

#combined phenotypes thr 0.05
ggplot_In3RMo_thr_0.05 <- ggplot(Enrichment_FULL_INV_3RMo_Rbound_thr0.05, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
  xlab("Chromosome Section") +
  ylab("nTOP") +
  geom_hline(yintercept=0.05, linetype = "dashed", color ="red") +
  facet_wrap(~thr) + 
  theme(axis.text.x=element_text(size=9)) +
  ggtitle("3RMo Enrichment, thr = 0.05")


'thr 0.01~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#null hypothesis that inv would have 0.001
Enrichment_FULL_INV_3RMo_thr0.05 <- foreach(i=1:length(VectorOf_3RMo_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  
  data <- fread(VectorOf_3RMo_Pheno[i])
  
  
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.01)), thr = c(0.01), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_3RMo_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
Enrichment_FULL_INV_3RMo_Rbound_thr0.01 <- rbindlist(Enrichment_FULL_INV_3RMo_thr0.01)

#combined phenotypes thr 0.01 
ggplot_In3RMo_thr_0.01 <- ggplot(Enrichment_FULL_INV_3RMo_Rbound_thr0.01, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
  xlab("inv") +
  ylab("nTOP") +
  geom_hline(yintercept=0.01, linetype = "dashed", color ="red") +
  facet_wrap(~thr) +
  theme(axis.text.x=element_text(size=9)) +
  ggtitle("3RMo Enrichment, thr = 0.01")


'thr 0.001~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
Enrichment_FULL_INV_3RMo_thr0.001 <- foreach(i=1:length(VectorOf_3RMo_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  
  data <- fread(VectorOf_3RMo_Pheno[i])
  
  
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR,inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_3RMo_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
Enrichment_FULL_INV_3RMo_Rbound_thr0.001 <- rbindlist(Enrichment_FULL_INV_3RMo_thr0.001)

#combined phenotypes thr 0.001 
ggplot_In3RMo_thr_0.001 <- ggplot(Enrichment_FULL_INV_3RMo_Rbound_thr0.001, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar") +
  xlab("inv") +
  ylab("nTOP") +
  geom_hline(yintercept=0.001, linetype = "dashed", color ="red") +
  facet_wrap(~thr) +
  theme(axis.text.x=element_text(size=9)) +
  ggtitle("3RMo Enrichment, thr = 0.001")


'all thr graphs, combined'
library(gridExtra)
# Arrange the plots in a grid with a 1x3 layout
grid.arrange(ggplot_In3RMo_thr_0.05, ggplot_In3RMo_thr_0.01, ggplot_In3RMo_thr_0.001, ncol = 1) 
#+grid.text("Tag A", x = unit(0.05, "npc"), y = unit(0.95, "npc"), just = "left", gp = gpar(fontsize = 16, fontface = "bold"))

'for all pheno, separated'
#per individual phenotype
ggplot(Enrichment_FULL_INV_3RMo_Rbound, aes(x=inv, y=nTop, fill = inv)) +
  stat_summary(fun = "mean", geom = "bar") +
  xlab("inv") +
  ylab("nTOP") +
  facet_wrap(~Phenotype) +
  geom_hline(yintercept=0.0010, linetype = "dashed", color ="red") +
  ggtitle("All Phenotypes Inversion 3RMo Enrichment Compared to Null Hypothesis (nTop = 0.001)")

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#Loop to graph GWAS hits with CHR position buckets
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  'In2Lt GWAS hits, thr = 0.001'################################################
  
PhenoGather_2Lt <- foreach(i=1:length(VectorOf_2Lt_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  #setwd("/Users/alanbergland/Daria Pheno")
  
  #data = VectorOfPheno[i] #set data
  #data = as.data.table(read.delim(data)) #loads data
  data <- fread(VectorOf_2Lt_Pheno[i])
  
  #data[,pa:=p.adjust(PVAL)]
  #data = subset(data, PVAL < 0.05)
  #data = subset(data, CHR == "2L"| CHR =="3R")
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, pos_bin, inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_2Lt_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
gather_2Lt = rbindlist(PhenoGather_2Lt)

#Chill coma pheno
gather_2Lt_ChillComaPheno = subset(gather_2Lt, Phenotype == "ChillComaRecoveryTime_standard_male.nogrms.original.txt" | Phenotype == "ChillComaRecoveryTime_standard_female.nogrms.original.txt")
  ggplot_gather_2Lt_ChillComaPheno <- ggplot(gather_2Lt_ChillComaPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Chill Coma Recovery Pheno ~~~ In(2L)t")
#Climbing pheno
gather_2Lt_ClimbingPheno = subset(gather_2Lt, Phenotype == "ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt" | Phenotype == "ClimbingHeight_standard_M.nogrms.jan.txt")
  ggplot_gather_2Lt_ClimbingPheno <- ggplot(gather_2Lt_ClimbingPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Climbing Pheno ~~~ In(2L)t")
#DaySleep Pheno
gather_2Lt_DaySleepPheno = subset(gather_2Lt, Phenotype == "DaySleep_standard_Female.nogrms.original.txt" | Phenotype == "DayBoutNumber_standard_Female.nogrms.original.txt")
  ggplot_gather_2Lt_DaySleepPheno <- ggplot(gather_2Lt_DaySleepPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Day-Sleep Pheno ~~~ In(2L)t")
#FoodIntake pheno
gather_2Lt_FoodIntake = subset(gather_2Lt, Phenotype == "FoodIntake_standard_female.nogrms.original.txt")
  ggplot_gather_2Lt_FoodIntake <- ggplot(gather_2Lt_FoodIntake, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("FoodIntake Pheno ~~~ In(2L)t")
#Glycerol + glucose diet pheno
gather_2Lt_GlycerolSet = subset(gather_2Lt, Phenotype == "FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt" | Phenotype =="FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt"| Phenotype == "FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt")
  ggplot_gather_2Lt_GlycerolSet <- ggplot(gather_2Lt_GlycerolSet, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Glycerol Set ~~~ In(2L)t")
#Ethanol exposures pheno
gather_2Lt_EthanolExposures = subset(gather_2Lt, Phenotype == "MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt" | Phenotype =="MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt")
  ggplot_gather_2Lt_EthanolExposures <- ggplot(gather_2Lt_EthanolExposures, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("EthanolExpores Set (M/F) ~~~ In(2L)t")
#mCT pheno
gather_2Lt_mCT = subset(gather_2Lt, Phenotype == "MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt" | Phenotype =="MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt" | Phenotype == "MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt")
  ggplot_gather_2Lt_mCT <- ggplot(gather_2Lt_mCT, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("mCT set - 29, 26, 23 degree C  ~~~ In(2L)t")
#Soluble protein level pheno
gather_2Lt_SolubleProteinLevels = subset(gather_2Lt, Phenotype == "SolubleProteinLevels_PooledDiets_male.nogrms.original.txt" | Phenotype =="SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt")
  ggplot_gather_2Lt_SolubleProteinLevels <- ggplot(gather_2Lt_SolubleProteinLevels, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("SolubleProteinLevels  ~~~ In(2L)t")
#Weight w glucose levels pheno
gather_2Lt_WeightwGlucoseLevles = subset(gather_2Lt, Phenotype == "Weight_PooledDiets_male.nogrms.original.txt" | Phenotype =="Weight_LowGlucoseDiet_male.nogrms.original.txt" | Phenotype == "Weight_HighGlucoseDiet_male.nogrms.original.txt")
  ggplot_gather_2Lt_WeightwGlucoseLevles <- ggplot(gather_2Lt_WeightwGlucoseLevles, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("WeightwGlucoseLevels  ~~~ In(2L)t")


'In3RMo GWAS hits, thr = 0.001' ################################################
PhenoGather_3RMo <- foreach(i=1:length(VectorOf_3RMo_Pheno)) %dopar% {
  #i = 1
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  #setwd("/Users/alanbergland/Daria Pheno")
  
  #data = VectorOfPheno[i] #set data
  #data = as.data.table(read.delim(data)) #loads data
  data <- fread(VectorOf_3RMo_Pheno[i])
  
  #data[,pa:=p.adjust(PVAL)]
  #data = subset(data, PVAL < 0.05)
  #data = subset(data, CHR == "2L"| CHR =="3R")
  data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
  data[,pos_bin:=round(POS, -4)]
  
  data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  #data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, PVAL, pa, rnp, pos_bin, inv)]
  data2 <- data[,list(nTop = c(mean(rnp < 0.001)), thr = c(0.001), n=length(rnp)), list(CHR, pos_bin, inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,Phenotype:= VectorOf_3RMo_Pheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
gather_3RMo = rbindlist(PhenoGather_3RMo)
    
#scent attraction pheno                                
gather_3RMo_ScentAttractionPheno = subset(gather_3RMo, Phenotype == "2-HeptanoneAttraction_standard_male.nogrms.original.txt" | Phenotype == "2-HeptanoneAttraction_standard_female.nogrms.original.txt" | Phenotype == "1-HexanolAttraction_standard_female.nogrms.original.txt" | Phenotype == "HexanalAttraction_standard_male.nogrms.original.txt" | Phenotype == "HexanalAttraction_standard_female.nogrms.original.txt")
  ggplot_gather_3RMo_ScentAttractionPheno <- ggplot(gather_3RMo_ScentAttractionPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Scent Attraction Pheno ~~~ In(3R)Mo")
#centroid pheno
gather_3RMo_CentroidPheno = subset(gather_3RMo, Phenotype == "CentroidSize_standard_F.nogrms.jan.txt" | Phenotype == "CentroidSize_standard_M.nogrms.jan.txt" | Phenotype == "WingCentroidSize_standard_female.nogrms.original.txt" | Phenotype == "WingCentroidSize_standard_male.nogrms.original.txt")
  ggplot_gather_3RMo_CentroidPheno <- ggplot(gather_3RMo_CentroidPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Centroid Pheno ~~~ In(3R)Mo")
#internocular distance pheno
gather_3RMo_InternocularDistancePheno = subset(gather_3RMo, Phenotype == "InternocularDistance_standard_F.nogrms.jan.txt" | Phenotype == "InternocularDistance_standard_M.nogrms.jan.txt")
  ggplot_gather_3RMo_InternocularDistancePheno <- ggplot(gather_3RMo_InternocularDistancePheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Internocular Distance Pheno ~~~ In(3R)Mo")
#triglyceride levels pheno
gather_3RMo_TriglycerideLevelPheno = subset(gather_3RMo, Phenotype == "TriglycerideLevels_PooledDiets_male.nogrms.original.txt" | Phenotype == "TriglycerideLevels_LowGlucoseDiet_male.nogrms.original.txt")
  ggplot_gather_3RMo_TriglycerideLevelPheno <- ggplot(gather_3RMo_TriglycerideLevelPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Triglyceride Level Pheno ~~~ In(3R)Mo")=
#weight/glucose pheno
gather_3RMo_WeightGlucosePheno = subset(gather_3RMo, Phenotype == "Weight_PooledDiets_male.nogrms.original.txt" | Phenotype == "Weight_HighGlucoseDiet_male.nogrms.original.txt" | Phenotype == "Weight_LowGlucoseDiet_male.nogrms.original.txt")
  ggplot_gather_3RMo_WeightGlucosePheno <- ggplot(gather_3RMo_WeightGlucosePheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Weight Glucose Pheno ~~~ In(3R)Mo")
#misc
gather_3RMo_MiscPheno = subset(gather_3RMo, Phenotype == "IntensityOfResponse_standard_F.nogrms.original.txt" | Phenotype == "ThoraxLength_standard_female.nogrms.original.txt" | Phenotype == "PupalCaseLength_Standard_Mix.nogrms.jan.txt" | Phenotype == "NegativeGeotaxis_MSB_Treatment_female.nogrms.original.txt" | Phenotype == "LifeTimeFecundity_LowYeast_female.nogrms.original.txt")
  ggplot_gather_3RMo_MiscPheno <- ggplot(gather_3RMo_MiscPheno, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phenotype~CHR) + ggtitle("Misc sPheno ~~~ In(3R)Mo")


      
#Testing whether overlaps between phenotypes are greater than expected 
VectorOfPheno <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "DaySleep_standard_Female.nogrms.original.txt", "Weight_PooledDiets_male.nogrms.original.txt") #add the rest of the file names here 
library(dplyr)
example_recombinations <- as.data.table(as.matrix(combn(VectorOfPheno, 2)))

'Example enrichment test'
setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
PhenoGather_Enrichment <- foreach(comb=example_recombinations, .combine='rbind') %do% {
    i <- comb[1]
    data1 <- fread(i)
    data1[,rnp_1:=rank(PVAL)/(length(PVAL)+1)]
    data1[,pos_bin:=round(POS, -4)]
    
    data1[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
    data1[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
    data1[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
    data1[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
    
    j = comb[2]
    data2 <- fread(j)
    data2[,rnp_2:=rank(PVAL)/(length(PVAL)+1)]
    data2[,pos_bin:=round(POS, -4)]
    
    setkey(data1, "CHR", "POS")
    setkey(data2, "CHR", "POS")
    data <- merge(data1, data2, by = c("CHR", "POS"))
    
    data2 <- data[,list(nTop = mean(rnp_1 < 0.05 & rnp_2 < .05), thr = c(0.05^2), n=length(rnp_1), eq=mean(rnp_1==rnp_2)),
                  list(CHR, pos_bin=pos_bin.x, inv = inv)]
    data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
    data2[,PhenotypeOne:= i] #can use a split to take off file extensions from names
    data2[,PhenotypeTwo:= j]
    return(data2)
}

ggplot(data=PhenoGather_Enrichment, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(.~CHR)


'Enrichment Test - In2Lt'
In2lt_recombinations <- as.data.table(as.matrix(combn(VectorOf_2Lt_Pheno, 2)))
In2Lt_Enrichment <- foreach(comb=In2lt_recombinations, .combine='rbind') %do% {
  
  i <- comb[1]
  data1 <- fread(i)
  data1[,rnp_1:=rank(PVAL)/(length(PVAL)+1)]
  data1[,pos_bin:=round(POS, -4)]
  
  data1[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data1[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data1[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data1[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  j = comb[2]
  data2 <- fread(j)
  data2[,rnp_2:=rank(PVAL)/(length(PVAL)+1)]
  data2[,pos_bin:=round(POS, -4)]
  
  setkey(data1, "CHR", "POS")
  setkey(data2, "CHR", "POS")
  
  data <- merge(data1, data2, by = c("CHR", "POS"))
  
  data2 <- data[,list(nTop = mean(rnp_1 < 0.05 & rnp_2 < .05), thr = c(0.05^2), n=length(rnp_1), eq=mean(rnp_1==rnp_2)),
                list(CHR, pos_bin=pos_bin.x, inv = inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,PhenotypeOne:= i] #can use a split to take off file extensions from names
  data2[,PhenotypeTwo:= j]
  return(data2)
}

In2Lt_Enrichment_plot <- ggplot(data=In2Lt_Enrichment, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(.~CHR) + ggtitle("In2Lt_Enrichment_plot")
write.csv(In2Lt_Enrichment, "/Users/dgg/Desktop\\Recomb_Bucket_Pvalues_In2Lt.csv", row.names=FALSE)

setwd("/Users/dgg/Desktop/Thesis/Recomb Files")
In2Lt_Enrichment <- fread("Recomb_Bucket_Pvalues_In2Lt.csv")


'Enrichment Test - In3RMo'
In3RMo_recombinations <- as.data.table(as.matrix(combn(VectorOf_3RMo_Pheno, 2)))
In3RMo_Enrichment <- foreach(comb=In3RMo_recombinations, .combine='rbind') %do% {
  
  i <- comb[1]
  data1 <- fread(i)
  data1[,rnp_1:=rank(PVAL)/(length(PVAL)+1)]
  data1[,pos_bin:=round(POS, -4)]
  
  data1[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data1[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data1[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data1[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]
  
  j = comb[2]
  data2 <- fread(j)
  data2[,rnp_2:=rank(PVAL)/(length(PVAL)+1)]
  data2[,pos_bin:=round(POS, -4)]
  
  setkey(data1, "CHR", "POS")
  setkey(data2, "CHR", "POS")
  
  data <- merge(data1, data2, by = c("CHR", "POS"))
  
  data2 <- data[,list(nTop = mean(rnp_1 < 0.05 & rnp_2 < .05), thr = c(0.05^2), n=length(rnp_1), eq=mean(rnp_1==rnp_2)),
                list(CHR, pos_bin=pos_bin.x, inv = inv)]
  
  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
  
  data2[,PhenotypeOne:= i] #can use a split to take off file extensions from names
  data2[,PhenotypeTwo:= j]
  return(data2)
}
In3RMo_Enrichment_plot <- ggplot(data=In3RMo_Enrichment, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(.~CHR) + ggtitle("In3RMo_Enrichment_plot")
write.csv(In3RMo_Enrichment, "/Users/dgg/Desktop\\Recomb_Bucket_Pvalues_In2Lt.csv", row.names=FALSE)

setwd("/Users/dgg/Desktop/Thesis/Recomb Files")
In3RMo_Enrichment <- fread("Recomb_Bucket_Pvalues_In3RMo.csv")


'~~~~~~~~~~~~~~~~~~~~~~~~~~~Background exploration into total significant densities and score~~~~~~~~~~~~~~~~~~~~~~~`'
  #background data
  In(2L)t --> position 2225744 (distal breakpoint = far from the centromere)
      chromosome 2L is --- to 2225987 (2225744 to 13154180, Adams Data
  In(3R)Mo --> position 17232639 (proximal breakpoint = close to the centromere)

  #splitting dataset by chromsomes 
    #In(2L)t
  CentroidSize_Standard_F_2Lt <- subset(CentroidSize_Standard_F, CHR == "2L")  #max = 22986299
    CentroidSize_Standard_F_2Lt_Wide <- subset(CentroidSize_Standard_F_2Lt, PVAL <=0.05)
  CentroidSize_Standard_M_2Lt <- subset(CentroidSize_Standard_M, CHR == "2L") 
    CentroidSize_Standard_M_2Lt_Wide <- subset(CentroidSize_Standard_M_2Lt, PVAL <=0.05)
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt <- subset(ClimbingHeight_ParaquatInsecticideExposure_M, CHR == "2L") 
    ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt, PVAL <=0.05)
  ClimbingHeight_standard_M_2Lt <- subset(ClimbingHeight_standard_M, CHR == "2L") 
    ClimbingHeight_standard_M_2Lt_Wide<- subset(ClimbingHeight_standard_M_2Lt, PVAL <=0.05)
  InternocularDistance_standard_F_2Lt <- subset(InternocularDistance_standard_F, CHR == "2L") 
    InternocularDistance_standard_F_2Lt_Wide<- subset(InternocularDistance_standard_F_2Lt, PVAL <=0.05)
  InternocularDistance_standard_M_2Lt <- subset(InternocularDistance_standard_M, CHR == "2L") 
    InternocularDistance_standard_M_2Lt_Wide<- subset(InternocularDistance_standard_M_2Lt, PVAL <=0.05)
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M, CHR == "2L") 
    MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide<- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt, PVAL <=0.05)
  MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M, CHR == "2L") 
    MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide<- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt, PVAL <=0.05)
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M, CHR == "2L") 
    MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide<- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt, PVAL <=0.05)
  PupalCaseLength_Standard_Mix_2Lt <- subset(PupalCaseLength_Standard_Mix, CHR == "2L") 
    PupalCaseLength_Standard_Mix_2Lt_Wide<- subset(PupalCaseLength_Standard_Mix_2Lt, PVAL <=0.05)
  
    #In(3R)Mo
  CentroidSize_Standard_F_3RMo <- subset(CentroidSize_Standard_F, CHR == "3R")  #max = 22986299
    CentroidSize_Standard_F_3RMo_Wide <- subset(CentroidSize_Standard_F_3RMo, PVAL <=0.05)
  CentroidSize_Standard_M_3RMo <- subset(CentroidSize_Standard_M, CHR == "3R") 
    CentroidSize_Standard_M_3RMo_Wide <- subset(CentroidSize_Standard_M_3RMo, PVAL <=0.05)
  ClimbingHeight_ParaquatInsecticideExposure_M_3RMo <- subset(ClimbingHeight_ParaquatInsecticideExposure_M, CHR == "3R") 
    ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo, PVAL <=0.05)
  ClimbingHeight_standard_M_3RMo <- subset(ClimbingHeight_standard_M, CHR == "3R") 
    ClimbingHeight_standard_M_3RMo_Wide<- subset(ClimbingHeight_standard_M_3RMo, PVAL <=0.05)
  InternocularDistance_standard_F_3RMo <- subset(InternocularDistance_standard_F, CHR == "3R") 
    InternocularDistance_standard_F_3RMo_Wide<- subset(InternocularDistance_standard_F_3RMo, PVAL <=0.05)
  InternocularDistance_standard_M_3RMo <- subset(InternocularDistance_standard_M, CHR == "3R") 
    InternocularDistance_standard_M_3RMo_Wide<- subset(InternocularDistance_standard_M_3RMo, PVAL <=0.05)
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M, CHR == "3R") 
    MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide<- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo, PVAL <=0.05)
  MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M, CHR == "3R") 
    MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide<- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo, PVAL <=0.05)
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M, CHR == "3R") 
    MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide<- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo, PVAL <=0.05)
  PupalCaseLength_Standard_Mix_3RMo <- subset(PupalCaseLength_Standard_Mix, CHR == "3R") 
    PupalCaseLength_Standard_Mix_3RMo_Wide<- subset(PupalCaseLength_Standard_Mix_3RMo, PVAL <=0.05)
  
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~~~~~P-value density graphs, all pheno, In(2L)t and In(3R)Mo~~~~~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
AHHHH <- ggplot() +
  geom_density(CentroidSize_Standard_F_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "purple") +
  geom_vline(xintercept = 2205744, color="black") +
  geom_vline(xintercept = 13174180, color="black") + 
  ggtitle("P<0.5 In(2Lt), All Pheno")

'experiment <- CentroidSize_Standard_F_2Lt_Wide$POS
experiment <- as.data.frame(experiment)
name <- c('POS')
colnames(experiment) <- name
AHHHH2 <- ggplot() +
  geom_density(experiment, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "purple") +
  geom_vline(xintercept = 2205744, color="black") +
  geom_vline(xintercept = 13174180, color="black") + 
  ggtitle("P<0.5 In(2Lt), All Pheno")'

 AllPheno_2Lt_Sig <- ggplot() +
    geom_density(CentroidSize_Standard_F_2Lt_Wide, mapping = aes(x=POS, y = (..count../sum(..count..))), color = "purple") +
    geom_density(CentroidSize_Standard_M_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "purple", linetype ="dashed") +
    geom_density(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "black") +
    geom_density(ClimbingHeight_standard_M_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "black", linetype ="dashed") +
    geom_density(InternocularDistance_standard_F_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "blue") +
    geom_density(InternocularDistance_standard_M_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "blue", linetype ="dashed") +
    geom_density(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "red") +
    geom_density(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide , mapping = aes(x=POS, y = ..count../sum(..count..)), color = "red", linetype = "dashed") +
    geom_density( MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide , mapping = aes(x=POS, y = ..count../sum(..count..)), color = "red", linetype = "dotted") +
    geom_density(PupalCaseLength_Standard_Mix_2Lt_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "black") +
    geom_vline(xintercept = 2205744, color="black") +
    geom_vline(xintercept = 13174180, color="black") + 
    labs(y = "Density, 1 = 100%", x = "Chromosome Position") +
    ggtitle("P<0.05 In(2Lt), All Pheno") +
    theme(plot.title = element_text(hjust = 0.5))
  
  AllPheno_3RMo_Sig <- ggplot() +
    geom_density(CentroidSize_Standard_F_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "purple") +
    geom_density(CentroidSize_Standard_M_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "purple", linetype ="dashed") +
    geom_density(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "black") +
    geom_density(ClimbingHeight_standard_M_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "black", linetype ="dashed") +
    geom_density(InternocularDistance_standard_F_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "blue") +
    geom_density(InternocularDistance_standard_M_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "blue", linetype ="dashed") +
    geom_density(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "red") +
    geom_density(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide , mapping = aes(x=POS, y = ..count../sum(..count..)), color = "red", linetype = "dashed") +
    geom_density( MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide , mapping = aes(x=POS, y = ..count../sum(..count..)), color = "red", linetype = "dotted") +
    geom_density(PupalCaseLength_Standard_Mix_3RMo_Wide, mapping = aes(x=POS, y = ..count../sum(..count..)), color = "black") +
    geom_vline(xintercept = 21406917, color="black") +
    geom_vline(xintercept = 29031297, color="black") + 
    geom_vline(xintercept = 16432209, color="yellow") +
    geom_vline(xintercept = 24744010, color="yellow") + 
    labs(y = "Density, 1 = 100%", x = "Chromosome Position") +
    ggtitle("P<0.05 In(3RMo), All Pheno") +
    theme(plot.title = element_text(hjust = 0.5))
  
  #individual graphs - left INV border - In(2L)t
  CentroidSize_Standard_F_2Lt_Trimmed <- subset(CentroidSize_Standard_F_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  CentroidSize_Standard_F_2Lt_LeftBorder <- ggplot(CentroidSize_Standard_F_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), CentroidSize_Standard_F")

  CentroidSize_Standard_M_2Lt_Trimmed <- subset(CentroidSize_Standard_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  CentroidSize_Standard_M_2Lt_Graph <- ggplot(CentroidSize_Standard_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), CentroidSize_Standard_M")
  
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Trimmed <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Graph <- ggplot(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), ClimbingHeight_ParaquatInsecticideExposure_M")
  
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt, PVAL <=0.05)
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_ALlInv <- 
    ggplot(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide, aes (x=POS)) +
    geom_density() +  
    geom_vline(xintercept = 2205744, color="purple") +
    geom_vline(xintercept = 13174180, color="purple") + 
    ggtitle("P < 0.5 In(2Lt), ClimbingHeight_ParaquatInsecticideExposure_M")
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_ALlInv
  
  ClimbingHeight_standard_M_2Lt_Trimmed <- subset(ClimbingHeight_standard_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  ClimbingHeight_standard_M_2Lt_Graph <- ggplot(ClimbingHeight_standard_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), ClimbingHeight_standard_M")
  
  ClimbingHeight_standard_M_2Lt_Wide <- subset(ClimbingHeight_standard_M_2Lt, PVAL <=0.05)
  ClimbingHeight_ClimbingHeight_standard_M_2Lt_ALlInv <- 
    ggplot(ClimbingHeight_standard_M_2Lt_Wide, aes (x=POS)) +
    geom_density() +  
    geom_vline(xintercept = 2205744, color="purple") +
    geom_vline(xintercept = 13174180, color="purple") + 
    ggtitle("P < 0.5 In(2Lt), ClimbingHeight_standard_M")
  ClimbingHeight_ClimbingHeight_standard_M_2Lt_ALlInv
  
  InternocularDistance_standard_F_2Lt_Trimmed <- subset(InternocularDistance_standard_F_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  InternocularDistance_standard_F_2Lt_Graph <- ggplot(InternocularDistance_standard_F_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), InternocularDistance_standard_F")
  
  InternocularDistance_standard_M_2Lt_Trimmed <- subset(InternocularDistance_standard_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  InternocularDistance_standard_M_2Lt_Graph <- ggplot(InternocularDistance_standard_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), InternocularDistance_standard_M")
  
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Trimmed <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Graph <- ggplot(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), mCT_TwentyNineDegreesC_M")
  
  MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Trimmed <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Graph <- ggplot(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), mCT_TwentySixDegreesC_M")
  
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Trimmed <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Graph <- ggplot(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), mCT_TwentyThreeDegreesC_M")
  
  PupalCaseLength_Standard_Mix_2Lt_Trimmed <- subset(PupalCaseLength_Standard_Mix_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  PupalCaseLength_Standard_Mix_2Lt_Graph <- ggplot(PupalCaseLength_Standard_Mix_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), PupalCaseLength_M")

  
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'~~~~~~Score & P-value density graphs, all pheno, In(2L)t and In(3R)Mo~~~~~~~~~~'
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
  
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  'In2Lt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  #prepping pos and neg score datsets
CentroidSize_Standard_F_2Lt_Wide_PosScore <- subset(CentroidSize_Standard_F_2Lt_Wide, SCORE > 0)
CentroidSize_Standard_F_2Lt_Wide_NegScore <- subset(CentroidSize_Standard_F_2Lt_Wide, SCORE < 0)

CentroidSize_Standard_M_2Lt_Wide_PosScore <- subset(CentroidSize_Standard_M_2Lt_Wide, SCORE > 0)
CentroidSize_Standard_M_2Lt_Wide_NegScore <- subset(CentroidSize_Standard_M_2Lt_Wide, SCORE < 0)

ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide_PosScore <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide, SCORE > 0)
ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide_NegScore <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide, SCORE < 0)

ClimbingHeight_standard_M_2Lt_Wide_PosScore <- subset(ClimbingHeight_standard_M_2Lt_Wide, SCORE > 0)
ClimbingHeight_standard_M_2Lt_Wide_NegScore <- subset(ClimbingHeight_standard_M_2Lt_Wide, SCORE < 0)

InternocularDistance_standard_F_2Lt_Wide_PosScore <- subset(InternocularDistance_standard_F_2Lt_Wide, SCORE > 0)
InternocularDistance_standard_F_2Lt_Wide_NegScore <- subset(InternocularDistance_standard_F_2Lt_Wide, SCORE < 0)

InternocularDistance_standard_M_2Lt_Wide_PosScore <- subset(InternocularDistance_standard_M_2Lt_Wide, SCORE > 0)
InternocularDistance_standard_M_2Lt_Wide_NegScore <- subset(InternocularDistance_standard_M_2Lt_Wide, SCORE < 0)

MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide_PosScore <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide, SCORE > 0)
MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide_NegScore <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide, SCORE < 0)

MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide_PosScore <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide, SCORE > 0)
MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide_NegScore <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide, SCORE < 0)

MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide_PosScore <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide, SCORE > 0)
MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide_NegScore <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide, SCORE < 0)

PupalCaseLength_Standard_Mix_2Lt_Wide_PosScore <- subset(PupalCaseLength_Standard_Mix_2Lt_Wide, SCORE > 0)
PupalCaseLength_Standard_Mix_2Lt_Wide_NegScore <- subset(PupalCaseLength_Standard_Mix_2Lt_Wide, SCORE < 0)

  'Score Graphs~~~~~~~In2Lt~~~~~~~~~~~~~~~~~~~~~'
CentroidSize_Standard_F_2Lt_Wide_Score_Graph <- ggplot() +
  geom_density(CentroidSize_Standard_F_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(CentroidSize_Standard_F_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(CentroidSize_Standard_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(CentroidSize_Standard_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(ClimbingHeight_ParaquatInsecticideExposure_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(ClimbingHeight_standard_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(ClimbingHeight_standard_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(InternocularDistance_standard_F_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(InternocularDistance_standard_F_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(InternocularDistance_standard_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(InternocularDistance_standard_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(PupalCaseLength_Standard_Mix_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(PupalCaseLength_Standard_Mix_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_vline(xintercept = 2205744, color="black") +
  geom_vline(xintercept = 13174180, color="black") + 
  labs(y = "Density", x = "Chromosome Position") +
  ggtitle("P<0.05 In(2Lt), All Pheno, score") +
  theme(plot.title = element_text(hjust = 0.5))

CentroidSize_Standard_F_2Lt_Wide_Score_Graph <- ggplot() +
  geom_density(CentroidSize_Standard_F_2Lt_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(CentroidSize_Standard_F_2Lt_Wide_NegScore, mapping = aes(x=POS), color = "red") + 
  geom_density(CentroidSize_Standard_F_2Lt_Wide, mapping = aes(x=POS), color = "black") +
  geom_vline(xintercept = 2205744, color="black") +
  geom_vline(xintercept = 13174180, color="black") + 
  ggtitle("P<0.05 In(2Lt), All Pheno, score") +
  theme(plot.title = element_text(hjust = 0.5))

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
'In3RMo~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#prepping pos and neg score datsets
CentroidSize_Standard_F_3RMo_Wide_PosScore <- subset(CentroidSize_Standard_F_3RMo_Wide, SCORE > 0)
CentroidSize_Standard_F_3RMo_Wide_NegScore <- subset(CentroidSize_Standard_F_3RMo_Wide, SCORE < 0)

CentroidSize_Standard_M_3RMo_Wide_PosScore <- subset(CentroidSize_Standard_M_3RMo_Wide, SCORE > 0)
CentroidSize_Standard_M_3RMo_Wide_NegScore <- subset(CentroidSize_Standard_M_3RMo_Wide, SCORE < 0)

ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide_PosScore <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide, SCORE > 0)
ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide_NegScore <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide, SCORE < 0)

ClimbingHeight_standard_M_3RMo_Wide_PosScore <- subset(ClimbingHeight_standard_M_3RMo_Wide, SCORE > 0)
ClimbingHeight_standard_M_3RMo_Wide_NegScore <- subset(ClimbingHeight_standard_M_3RMo_Wide, SCORE < 0)

InternocularDistance_standard_F_3RMo_Wide_PosScore <- subset(InternocularDistance_standard_F_3RMo_Wide, SCORE > 0)
InternocularDistance_standard_F_3RMo_Wide_NegScore <- subset(InternocularDistance_standard_F_3RMo_Wide, SCORE < 0)

InternocularDistance_standard_M_3RMo_Wide_PosScore <- subset(InternocularDistance_standard_M_3RMo_Wide, SCORE > 0)
InternocularDistance_standard_M_3RMo_Wide_NegScore <- subset(InternocularDistance_standard_M_3RMo_Wide, SCORE < 0)

MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide_PosScore <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide, SCORE > 0)
MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide_NegScore <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide, SCORE < 0)

MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide_PosScore <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide, SCORE > 0)
MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide_NegScore <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide, SCORE < 0)

MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide_PosScore <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide, SCORE > 0)
MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide_NegScore <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide, SCORE < 0)

PupalCaseLength_Standard_Mix_3RMo_Wide_PosScore <- subset(PupalCaseLength_Standard_Mix_3RMo_Wide, SCORE > 0)
PupalCaseLength_Standard_Mix_3RMo_Wide_NegScore <- subset(PupalCaseLength_Standard_Mix_3RMo_Wide, SCORE < 0)

  'Score Graphs~~~~~~~In3RMo~~~~~~~~~~~~~~~~~~~~~'
CentroidSize_Standard_F_3RMo_Wide_Score_Graph <- ggplot() +
  geom_density(CentroidSize_Standard_F_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(CentroidSize_Standard_F_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(CentroidSize_Standard_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(CentroidSize_Standard_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(ClimbingHeight_standard_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(ClimbingHeight_standard_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(InternocularDistance_standard_F_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(InternocularDistance_standard_F_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(InternocularDistance_standard_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(InternocularDistance_standard_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_density(PupalCaseLength_Standard_Mix_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(PupalCaseLength_Standard_Mix_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") +
  geom_vline(xintercept = 21406917, color="black") +
  geom_vline(xintercept = 29031297, color="black") + 
  geom_vline(xintercept = 16432209, color="yellow") +
  geom_vline(xintercept = 24744010, color="yellow") + 
  labs(y = "Density", x = "Chromosome Position") +
  ggtitle("P<0.05 In(3RMo), All Pheno") +
  theme(plot.title = element_text(hjust = 0.5))

CentroidSize_Standard_F_3RMo_Wide_Score_Graph <- ggplot() +
  geom_density(CentroidSize_Standard_F_3RMo_Wide_PosScore, mapping = aes(x=POS), color = "blue") +   
  geom_density(CentroidSize_Standard_F_3RMo_Wide_NegScore, mapping = aes(x=POS), color = "red") + 
  geom_density(CentroidSize_Standard_F_3RMo_Wide, mapping = aes(x=POS), color = "black") +
  geom_vline(xintercept = 21406917, color="black") +
  geom_vline(xintercept = 29031297, color="black") + 
  geom_vline(xintercept = 16432209, color="yellow") +
  geom_vline(xintercept = 24744010, color="yellow") + 
  labs(y = "Density by Percent", x = "Chromosome Position") +
  ggtitle("P<0.05 In(3RMo), All Pheno") +
  theme(plot.title = element_text(hjust = 0.5))

