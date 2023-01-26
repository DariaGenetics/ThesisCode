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

  #load in data
  setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  CentroidSize_Standard_F <- read.delim("CentroidSize_standard_F.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
  CentroidSize_Standard_M <- read.delim("CentroidSize_standard_M.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
  ClimbingHeight_ParaquatInsecticideExposure_M <- read.delim("ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
  ClimbingHeight_standard_M <- read.delim("ClimbingHeight_standard_M.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
  InternocularDistance_standard_F <- read.delim("InternocularDistance_standard_F.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
  InternocularDistance_standard_M <- read.delim("InternocularDistance_standard_M.nogrms.jan.txt",
                                        stringsAsFactor = FALSE)
  MinimumSurvivableTemperature_TwentyNineDegreesC_M <- read.delim("MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt",
                                                stringsAsFactor = FALSE)
  MinimumSurvivableTemperature_TwentySixDegreesC_M <- read.delim("MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt",
                                                                  stringsAsFactor = FALSE)
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M <- read.delim("MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt",
                                                                  stringsAsFactor = FALSE)
  PupalCaseLength_Standard_Mix <- read.delim("PupalCaseLength_Standard_Mix.nogrms.jan.txt",
                                                                  stringsAsFactor = FALSE)
  #background data
  In(2L)t --> position 2225744 (distal breakpoint = far from the centromere)
      chromosome 2L is --- to 2225987
  In(3R)Mo --> position 17232639 (proximal breakpoint = close to the centromere)
  
  #seeing how many data entries correspond to positions of note? 
  which(CentroidSize_Standard_F$POS == 2225744) #0
  which(CentroidSize_Standard_M$POS == 2225744) #0
  which(ClimbingHeight_ParaquatInsecticideExposure_M$POS == 2225744) #0
  which(ClimbingHeight_standard_M$POS == 2225744) #0
  which(InternocularDistance_standard_F$POS == 2225744) #0
  which(InternocularDistance_standard_M$POS == 2225744) #0
  which(MinimumSurvivableTemperature_TwentyNineDegreesC_M$POS == 2225744) #0
  which(MinimumSurvivableTemperature_TwentySixDegreesC_M$POS == 2225744) #0
  which(MinimumSurvivableTemperature_TwentyThreeDegreesC_M$POS == 2225744) #0
  which(PupalCaseLength_Standard_Mix$POS == 2225744) #0
  
  
  #splitting dataset by chromsomes 
    #In(2L)t
  CentroidSize_Standard_F_2Lt <- subset(CentroidSize_Standard_F, CHR == "2L") 
  CentroidSize_Standard_M_2Lt <- subset(CentroidSize_Standard_M, CHR == "2L") 
  ClimbingHeight_ParaquatInsecticideExposure_M_2Lt <- subset(ClimbingHeight_ParaquatInsecticideExposure_M, CHR == "2L") 
  ClimbingHeight_standard_M_2Lt <- subset(ClimbingHeight_standard_M, CHR == "2L") 
  InternocularDistance_standard_F_2Lt <- subset(InternocularDistance_standard_F, CHR == "2L") 
  InternocularDistance_standard_M_2Lt <- subset(InternocularDistance_standard_M, CHR == "2L") 
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_2Lt <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M, CHR == "2L") 
  MinimumSurvivableTemperature_TwentySixDegreesC_M_2Lt <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M, CHR == "2L") 
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_2Lt <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M, CHR == "2L") 
  PupalCaseLength_Standard_Mix_2Lt <- subset(PupalCaseLength_Standard_Mix, CHR == "2L") 
  
    #In(3R)Mo
  CentroidSize_Standard_F_3RMo <- subset(CentroidSize_Standard_F, CHR == "3R") 
  CentroidSize_Standard_M_3RMo <- subset(CentroidSize_Standard_M, CHR == "3R") 
  ClimbingHeight_ParaquatInsecticideExposure_M_3RMo <- subset(ClimbingHeight_ParaquatInsecticideExposure_M, CHR == "3R") 
  ClimbingHeight_standard_M_3RMo <- subset(ClimbingHeight_standard_M, CHR == "3R") 
  InternocularDistance_standard_F_3RMo <- subset(InternocularDistance_standard_F, CHR == "3R") 
  InternocularDistance_standard_M_3RMo <- subset(InternocularDistance_standard_M, CHR == "3R") 
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M, CHR == "3R") 
  MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M, CHR == "3R") 
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M, CHR == "3R") 
  PupalCaseLength_Standard_Mix_3RMo <- subset(PupalCaseLength_Standard_Mix, CHR == "3R") 
  
  #In(2L)t  PValue Graphs
  CentroidSize_Standard_F_2Lt_Trimmed <- subset(CentroidSize_Standard_F_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  CentroidSize_Standard_F_2Lt_Graph <- ggplot(CentroidSize_Standard_F_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
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
  
  ClimbingHeight_standard_M_2Lt_Trimmed <- subset(ClimbingHeight_standard_M_2Lt, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  ClimbingHeight_standard_M_2Lt_Graph <- ggplot(ClimbingHeight_standard_M_2Lt_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(2Lt), ClimbingHeight_standard_M")
  
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
  
  #In(3R)Mo  PValue Graphs
  
  CentroidSize_Standard_F_3RMo_Trimmed <- subset(CentroidSize_Standard_F_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  CentroidSize_Standard_F_3RMo_Graph <- ggplot(CentroidSize_Standard_F_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), CentroidSize_Standard_F")
  
  
  CentroidSize_Standard_M_3RMo_Trimmed <- subset(CentroidSize_Standard_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  CentroidSize_Standard_M_3RMo_Graph <- ggplot(CentroidSize_Standard_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), CentroidSize_Standard_M")
  
  ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Trimmed <- subset(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Graph <- ggplot(ClimbingHeight_ParaquatInsecticideExposure_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), ClimbingHeight_ParaquatInsecticideExposure_M")
  
  ClimbingHeight_standard_M_3RMo_Trimmed <- subset(ClimbingHeight_standard_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  ClimbingHeight_standard_M_3RMo_Graph <- ggplot(ClimbingHeight_standard_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), ClimbingHeight_standard_M")
  
  InternocularDistance_standard_F_3RMo_Trimmed <- subset(InternocularDistance_standard_F_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  InternocularDistance_standard_F_3RMo_Graph <- ggplot(InternocularDistance_standard_F_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), InternocularDistance_standard_F")
  
  InternocularDistance_standard_M_3RMo_Trimmed <- subset(InternocularDistance_standard_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  InternocularDistance_standard_M_3RMo_Graph <- ggplot(InternocularDistance_standard_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), InternocularDistance_standard_M")
  
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Trimmed <- subset(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Graph <- ggplot(MinimumSurvivableTemperature_TwentyNineDegreesC_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), mCT_TwentyNineDegreesC_M")
  
  MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Trimmed <- subset(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Graph <- ggplot(MinimumSurvivableTemperature_TwentySixDegreesC_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), mCT_TwentySixDegreesC_M")
  
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Trimmed <- subset(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Graph <- ggplot(MinimumSurvivableTemperature_TwentyThreeDegreesC_M_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), mCT_TwentyThreeDegreesC_M")
  
  PupalCaseLength_Standard_Mix_3RMo_Trimmed <- subset(PupalCaseLength_Standard_Mix_3RMo, POS > 2215000 & POS < 2237000) # In(2L)t --> position 2225744, max pos = 22986299
  PupalCaseLength_Standard_Mix_3RMo_Graph <- ggplot(PupalCaseLength_Standard_Mix_3RMo_Trimmed, aes (x=POS, y=PVAL)) +
    geom_point() +  
    geom_vline(xintercept = 2225744, color="purple") +
    geom_hline(yintercept = 0.05, color="red") +
    ggtitle("PValues for Positions near In(3RMo), PupalCaseLength_M")
  
  
  
  CentroidSize_Standard_F[CentroidSize_Standard_F$POS>2225500 & CentroidSize_Standard_F$POS<2226000, ]
  CentroidSize_Standard_F[CentroidSize_Standard_F$POS>2225300 & CentroidSize_Standard_F$POS<2225830, ]
 