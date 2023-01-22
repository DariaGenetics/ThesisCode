#Note - project in 'ThesisCode' repository

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

'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~In2LtStatus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
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

'~~~PCA analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
### libraries

### load data
MeltedPrincipleComponents <- as.data.table(melt(pc_out,
                                                id.var = c("ral_id","In2Lt"),
                                                variable.name = "PC", value.name="ComponentValues"))

### iterate through pcs
A <- foreach(i=unique(MeltedPrincipleComponents$PC), .combine="rbind") %dopar% {
  message(i)
  #i = unique(MeltedPrincipleComponents$PC)[100]
  
  Xa <- lm(ComponentValues ~ In2Lt, data = MeltedPrincipleComponents[PC == i])
  
  foo <- boot(MeltedPrincipleComponents[PC == i],
              function(data, indices) summary(lm(ComponentValues~In2Lt, data[indices,]))$r.squared,
              R=1000)
  
  
  out <- data.table(i=i, PC=as.numeric(gsub("Dim.", "", i)),
                    lci=quantile(foo$t, c(0.025)),
                    uci=quantile(foo$t, c(0.975)),
                    r2=summary(Xa)$r.squared,
                    p=summary(aov(Xa))[[1]][1,5],
                    pve=as.numeric(as.data.table(pc$eig)[as.numeric(gsub("Dim.", "", i)),"percentage of variance", with=T]))
  
  return(out)
}
A[,pa:=p.adjust(p)]

r2_plot <-
  ggplot(data=A, aes(x=PC, y=r2)) +
  geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
  geom_point() +
  geom_point(data=A[p<.05], aes(x=PC, y=r2), color="red") +
  geom_point(data=A[pa<.1], aes(x=PC, y=r2), color="blue") +
  ggtitle("Principle Components - IN2LT") + 
  theme(plot.title = element_text(size=12, hjust = 0.5))


### which factors contribute to Dim.8?
loading <- as.data.table(pc$var$coord)
loading[,pheno:=row.names(pc$var$coord)]
ll <- melt(loading, id.var="pheno")

ll.rank <- ll[,list(quan=rank(abs(value))/(length(value)+1), pheno=pheno, rank=rank(value)), list(variable)]

setkey(ll, variable, pheno)
setkey(ll.rank, variable, pheno)
ll <- merge(ll, ll.rank)

'full blue - 8'
'full red - 4, 9, 112, 128, 137, 145, 179, 191, 197'

'PC8 - Blue'
phenoLoad8 <-
  ggplot(data=ll[variable=="Dim.8"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.8"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
  ggtitle("Top 10% Influencing Phenotypes: PC 8
      p<0.05 - Bonferroni Correction") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

'PC4 - Red'
phenoLoad4 <-
  ggplot(data=ll[variable=="Dim.4"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.4"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
  ggtitle("Top 10% Influencing Phenotypes: PC 4
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

'PC9 - Red'
phenoLoad9 <-
  ggplot(data=ll[variable=="Dim.9"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.9"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype") +
  ggtitle("Top 10% Influencing Phenotypes: PC 9
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

'PC112 - Red'
phenoLoad112 <-
  ggplot(data=ll[variable=="Dim.112"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.112"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

'PC128 - Red'
phenoLoad128 <-
  ggplot(data=ll[variable=="Dim.128"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.128"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

'PC137 - Red'
phenoLoad137 <-
  ggplot(data=ll[variable=="Dim.137"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.137"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

'PC145 - Red'
phenoLoad145 <-
  ggplot(data=ll[variable=="Dim.145"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.145"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

'PC179 - Red'
phenoLoad179 <-
  ggplot(data=ll[variable=="Dim.179"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.179"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

'PC191 - Red'
phenoLoad191 <-
  ggplot(data=ll[variable=="Dim.191"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.191"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

'PC197 - Red'
phenoLoad197 <-
  ggplot(data=ll[variable=="Dim.197"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.197"][quan>.9],
            aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")

### effect plot
pc.ag <- MeltedPrincipleComponents[, list(mu=mean(ComponentValues), sd=sd(ComponentValues), n=length(ComponentValues)),
                                   list(In2Lt, PC)]
pc.ag[,se:=sd/sqrt(n)]

effect_plot8 <- ggplot(data=pc.ag[PC=="Dim.8"], aes(x=In2Lt, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se)) +
  ggtitle("Effect Size: PC 8
      p<0.05 - Bonferroni Correction") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

effect_plot4 <- ggplot(data=pc.ag[PC=="Dim.4"], aes(x=In2Lt, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se)) +
  ggtitle("Effect Size: PC 4
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

effect_plot9 <- ggplot(data=pc.ag[PC=="Dim.9"], aes(x=In2Lt, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se)) +
  ggtitle("Effect Size: PC 9
      p<0.05 - ANOVA Test") + 
  theme(plot.title = element_text(size=8, hjust = 0.5))

### mega plot
layout <- "
  AAAA
  BBCC"

r2_plot + phenoLoad8 + effect_plot8 +
  plot_layout(design=layout) + plot_annotation(tag_levels ="A")

layout <- "
  AAAA
  BBCC
  DDEE
  FFGG"

r2_plot + phenoLoad8 + effect_plot8 + phenoLoad4 + effect_plot4 + phenoLoad9 + effect_plot9 +
  plot_layout(design=layout) + plot_annotation(tag_levels ="A")



   
   
   
   
   
   
   

