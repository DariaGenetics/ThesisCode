library(ggplot2)
library(foreach)
library(data.table)
library(dplyr)
library(stringr)
library(ggpubr)

'~~~~~~~~~~~~99 bootstraps + SD for shifting whole inv location~~~~~~~~~~~~~~~~~'


VectorOf_2Lt_Pheno_OneFive <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "ChillComaRecoveryTime_standard_female.nogrms.original.txt",
                        "ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt", "ClimbingHeight_standard_M.nogrms.jan.txt",
                        "DaySleep_standard_Female.nogrms.original.txt") 
VectorOf_2Lt_Pheno_SixTwelve <- c("FoodIntake_standard_female.nogrms.original.txt", "FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt", 
                        "FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt", "FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt", 
                        "LarvaeSurvival_40μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt", "LarvaeSurvival_20μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt", 
                        "LarvaeSurvival_10μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt")

VectorOf_2Lt_Pheno_ThirteenEighteen <- c("MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt", "MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt", 
                        "MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt", "MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt", 
                        "MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt", "SolubleProteinLevels_PooledDiets_male.nogrms.original.txt")
                        
                        
VectorOf_2Lt_Pheno_NineteenTwentyTwo <- c("SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt", "Weight_PooledDiets_male.nogrms.original.txt", 
                        "Weight_LowGlucoseDiet_male.nogrms.original.txt", "Weight_HighGlucoseDiet_male.nogrms.original.txt") 
        

# Create an empty list to hold the results
results <- list()

# Loop through each file in 'VectorOf_2Lt_PhenoTest'
for (j in seq_along(VectorOf_2Lt_Pheno_SixTwelve)) {
  # Loop 100 times
  for (i in 1:100) {
    # Set the working directory and read in the data
    #i = 1
    setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
    datafull <- fread(VectorOf_2Lt_Pheno_Trial[j])
    data <- subset(datafull, CHR == "2L")
    
    # Perform the calculations and save the results to a variable
    data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
    data[,pos_bin:=round(POS, -4)]
    
    #left_inversion_boundary <- sample (1:16968732, 1)
    left_inversion_boundary <- if(i < 1.5) 2225744 else sample (1:16968732, 1) #CAUTION: if you use this loop again, adjust 1:16968732 to 1:12057863
    range_length <- 10928436
    
    data[CHR =="2L" & POS >= left_inversion_boundary & POS <= left_inversion_boundary + range_length, inv:="2Lt"]
    data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
    
    data2 <- data[, list(nTop = c(mean(rnp < 0.001)), SigNumber = c(sum(rnp <= 0.001)), thr = c(0.001), n = length(PVAL)), list(CHR,inv)]
    data2[,TotalFraction := SigNumber/n]                   
    
    data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]
    data2[,Phenotype:= VectorOf_2Lt_Pheno_SixTwelve[j]]
    data2[,Iteration:=i]
    
    # Append the results to the 'results' list
    results[[paste0("Phenotype", j ,"iteration", i)]] <- data2
  }
} 

# Combine the results into a single data frame
ResampledResults <- rbindlist(results)
write.csv(ResampledResults,"/Users/dgg/Enrichment6-12(p=0.001).csv")

setwd("/Users/dgg/Desktop") 
EnrichmentResamples <- fread("(p=0.05)EnrichmentBootstraps.csv")
    #(p=0.05)EnrichmentBootstraps.csv
    #(p=0.01)EnrichmentBootstraps.csv

ResampledResultsCulled <- subset(EnrichmentResamples, inv == "2Lt")

# Compute standard deviation by phenotype
ResampledResultsCulled_noIter1 <- ResampledResultsCulled %>% filter(Iteration > 1)
sd_by_phenotype <- ResampledResultsCulled_noIter1 %>%
  group_by(Phenotype) %>%
  summarize(mean_TotalFraction = mean(TotalFraction),
            sd = sd(TotalFraction))

# Merge standard deviation back into the main data frame
ResampledResultsCulled <- merge(ResampledResultsCulled, sd_by_phenotype, by="Phenotype", suffixes=c("", "_sd"))



p1 <- ggplot(ResampledResultsCulled, aes(x=inv, y=TotalFraction, color=factor(Iteration==1))) + 
  ggtitle("Contained Fraction of Total Significant SNPs In In2Lt (blue) and Resampled Comparisons (red) for Significance Threshold p < 0.01") +
  labs(y = "Total Fraction of SNPs", x = "Inversion 2Lt") +
  geom_errorbar(aes(ymin=mean_TotalFraction-sd, ymax=mean_TotalFraction+sd), width=0.1, color="black", position="identity") +
  geom_point(alpha=0.5, position = position_dodge(0)) +
  geom_point(aes(z=ifelse(Iteration==1, 1, 0)), position = position_dodge(0.2)) +
  scale_color_discrete(name="Iteration") +
  theme(legend.position="top") +
  facet_wrap(~Phenotype) 

NameVector = c("ChillComaRecoveryTime_M", "ChillComaRecoveryTime_F", 
               "ClimbingHeightInsecticide_M", "ClimbingHeight_M",
               "DaySleep_F",
               "FoodIntake_F", "GlycerolLevels_PooledDiets_M",
               "GlycerolLevels_LowGlucoseDiet_M","GlycerolLevels_HighGlucoseDiet_M",
               "LarvaeSurvival_40μg-mLCHL_Mix", "LarvaeSurvival_20μg-mLCHL_Mix",
               "LarvaeSurvival_10μg-mLCHL_Mix",
               "ElutionTime_EtOH_M", "ElutionTime_EtOH_F",
               "MinSurvivableTemp_29C_M", "MinSurvivableTemp_26C_M",
               "MinSurvivableTemp_23C_M", "SolutableProtein_PooledDiets_M",
               "SolutableProtein_LowGlucoseDiets_M","Weight_PooledDiets_M",
               "Weight_LowGlucoseDiets_M", "Weight_HighGlucoseDiets_M")

p2 + facet_grid(inv ~ Phenotype, labeller = labeller(Phenotype = NameVector))

'~~~~~~~~~~~~~~~~~~~Making more sophisticated with graph structure~~~~~~~~~~~~~~'

#scratchwork figuring out where the breakpoints are:
  #chromosome length
   #left boudary = 2225744, therefore 
    #LeftBreakpointLeftBorder = 2225744 - (1000000/2) = 1725744
    #LeftBreakpointRightBorder = LeftBreakpointLeftBorder + 1000000 
    #RightBreakpointLeftBorder = 13154180 - (1000000/2) = 12654180
    #RightBreakpointRightBorder = LeftBreakpointLeftBorder + 1000000 
    #InversionMiddle = 

VectorOf_2Lt_Pheno_Trial <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "ChillComaRecoveryTime_standard_female.nogrms.original.txt")
              
VectorOf_2Lt_Pheno_ThreeFive <- c("ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt", "ClimbingHeight_standard_M.nogrms.jan.txt", "DaySleep_standard_Female.nogrms.original.txt")                   

VectorOf_2Lt_Pheno_SixTwelve <- c("FoodIntake_standard_female.nogrms.original.txt", "FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt", 
                                  "FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt", "FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt", 
                                  "LarvaeSurvival_40μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt", "LarvaeSurvival_20μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt", 
                                  "LarvaeSurvival_10μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt")

VectorOf_2Lt_Pheno_ThirteenEighteen <- c("MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt", "MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt", 
                                         "MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt", "MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt", 
                                         "MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt", "SolubleProteinLevels_PooledDiets_male.nogrms.original.txt")
VectorOf_2Lt_Pheno_NineteenTwentyTwo <- c("SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt", "Weight_PooledDiets_male.nogrms.original.txt", 
                                          "Weight_LowGlucoseDiet_male.nogrms.original.txt", "Weight_HighGlucoseDiet_male.nogrms.original.txt") 

# Create an empty list to hold the results
ReshufflingResults <- list()

# Loop through each file in 'VectorOf_2Lt_PhenoTest'
for (j in seq_along(VectorOf_2Lt_Pheno_NineteenTwentyTwo)) {
  # Loop 100 times
  for (i in 1:101) {
    # Set the working directory and read in the data
    #i = 1
    setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
    datafull <- fread(VectorOf_2Lt_Pheno_NineteenTwentyTwo[j])
    data <- subset(datafull, CHR == "2L") #cuts on the amount of data, slides the breakpoints along this chromosome alone
    
    # Perform the calculations and save the results to a variable
    data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
    data[, Threshold := 0.001]
    #data[,pos_bin:=round(POS, -4)]
    
    #setting edges
    #if (LEFT EDGE LEFT BREAKPOINT TRUE) else if (LEFT EDGE RIGHT BREAKPOINT TRUE)) else if (SAMPLE)
    LeftBreakpointLeftBorder <- if(i < 1.5) 1725744 else sample (1:21986299, 1) #2225744 - (1000000/2) = inversion left edge - 1/2 Mb (megabase = 1 million basepairs), else length of chr 2L (22986299) - 1000000
    #RightBreakpointLeftBorder <- if(i < 1.5) 12654180 else sample (1:21986299, 1) #13154180 - (1000000/2) = inversion right edge - 1/2 Mb (megabase = 1 million basepairs), else length of chr 2L (22986299) - 1000000
    RightBreakpointLeftBorder <- if(i < 1.5) 12654180 else sample (1:21986299, 1) 
    BreakpointLength <- 1000000
    
    INVMiddleLeftBorder <- if(i<1.5) 2725744 else sample(1:13057863,1) #left bound actual: 2225744 + (1000000/2) = inversion left edge + 1/2 Mb (megabase = 1 million basepairs), else length of chr 2L (22986299) - 
    INVMiddleLength <- 9928436  # (13154180 - (1000000/2)) - (2225744 + (1000000/2)) = 9928436
   
    #FullINVLeftBorderForOutside <- if(i<1.5) 1725744 else sample(1:11057863, 1) #left breakpoint left edge to as far right as possible without full inversion + breakpoint length clipping out of chr. full length of chr = 22986299, width of breakpoint to breakpoint edge is 11928436
    #FullINVPLusBreakpoints <- 11928436
    
    data[POS >= LeftBreakpointLeftBorder & POS <= LeftBreakpointLeftBorder + BreakpointLength, LeftBreakpoint:="LeftBreakpoint"]
    data[POS >= RightBreakpointLeftBorder & POS <= RightBreakpointLeftBorder + BreakpointLength, RightBreakpoint:="RightBreakpoint"]
    data[POS >= INVMiddleLeftBorder & POS <= INVMiddleLeftBorder + INVMiddleLength, InversionMiddle:="InversionMiddle"]
    #data[POS <= FullINVLeftBorderForOutside & POS >= FullINVLeftBorderForOutside + FullINVPLusBreakpoints, OutsideInversion:="OutsideInversion"]
    
    #set(data, j="Threshold", value=Threshold)
    
    data2 <- data[, list(nTop_L_Breakpoint = c(mean(rnp < 0.001)), SigNumber_L_Breakpoint = c(sum(rnp <= 0.001)), n_L_Breakpoint = length(PVAL)), list(LeftBreakpoint)]
    data3 <- data[, list(nTop_R_Breakpoint = c(mean(rnp < 0.001)), SigNumber_R_Breakpoint = c(sum(rnp <= 0.001)), n_R_Breakpoint = length(PVAL)), list(RightBreakpoint)]
    data4 <- data[, list(nTop_Middle = c(mean(rnp < 0.001)), SigNumber_Middle = c(sum(rnp <= 0.001)), n_Middle = length(PVAL)), list(InversionMiddle)] 
    
    #data$nTopExp=mean(data[rnp<Threshold]$rnp,rm.na = T)
    data5 <- cbind(data2,data3,data4)
  
    data5[, Threshold := list(0.001)]
    data5[, Iteration := i]
    #data2[, p := pbinom(nTop * n, n, thr, lower.tail = F)] <---- is this necessary?
    data5[, Phenotype := VectorOf_2Lt_Pheno_NineteenTwentyTwo[j]]
    #data2[,PercentDifferenceFromThr:= ((nTop - thr)/((nTop + thr)/2))*100] 
    
    # Append the results to the 'results' list
    ReshufflingResults[[paste0("Phenotype", j ,"iteration", i)]] <- data5
  }
} 

# Combine the results into a single data frame
RBoundReshufflingResults <- rbindlist(ReshufflingResults)
write.csv(RBoundReshufflingResults,"/Users/dgg/TrialRunWith100_19-22.csv")
- 'TrialRunWith100_3-5.csv'

RBoundReshufflingResults <- fread("ThreeShuffles0.001.csv")


# Compute standard deviation by phenotype
ResampledResultsCulled_noIter1 <- RBoundReshufflingResults %>% filter(Iteration > 1) 

  LeftBreakpoint_sd_by_phenotype <- ResampledResultsCulled_noIter1 %>%
  group_by(Phenotype) %>%
  summarize(mean_nTop_LeftBreakpoint = mean(nTop_L_Breakpoint),
            sd_nTop_LeftBreakpoint = sd(nTop_L_Breakpoint))
  
  RightBreakpoint_sd_by_phenotype <- ResampledResultsCulled_noIter1 %>%
    group_by(Phenotype) %>%
    summarize(mean_nTop_RightBreakpoint = mean(nTop_R_Breakpoint),
              sd_nTop_RightBreakpoint = sd(nTop_R_Breakpoint))
  
  InvMiddle_sd_by_phenotype <- ResampledResultsCulled_noIter1 %>%
    group_by(Phenotype) %>%
    summarize(mean_nTop_InvMiddle = mean(nTop_Middle),
              sd_nTop_InvMiddle = sd(nTop_Middle))

# Merge standard deviation back into the main data frame
ResampledResultsCulled1 <- merge(RBoundReshufflingResults, LeftBreakpoint_sd_by_phenotype, by="Phenotype", suffixes=c("", "_sd"))
ResampledResultsCulled2 <- merge(ResampledResultsCulled1, RightBreakpoint_sd_by_phenotype, by="Phenotype", suffixes=c("", "_sd"))
ResampledResultsCulledFinal <- merge(ResampledResultsCulled2, InvMiddle_sd_by_phenotype, by="Phenotype", suffixes=c("", "_sd"))

#removes all "NA" rows
ResampledResultsCulledExperimental <- ResampledResultsCulledFinal[!(ResampledResultsCulledFinal$LeftBreakpoint=="NA")]
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'ChillComaRecoveryTime_standard_female.nogrms.original.txt'] <- 'ChillComaRecoveryTime_F'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'ChillComaRecoveryTime_standard_male.nogrms.original.txt'] <- 'ChillComaRecoveryTime_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt'] <- 'Climbing_InsecticideExposure_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'ClimbingHeight_standard_M.nogrms.jan.txt'] <- 'Climbing_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'DaySleep_standard_Female.nogrms.original.txt'] <- 'CircadianRhythm_F'
                 

'6-12'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'FoodIntake_standard_female.nogrms.original.txt'] <- 'FoodIntake_F'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt'] <- 'GlycerolLevels_PooledDiets_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt'] <- 'GLycerolLevels_LowGlucoseDiet_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt'] <- 'GLycerolLevels_HighGlucoseDiet_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'LarvaeSurvival_40μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt'] <- 'LarvaeSurvival_40μg-mLInsecticide_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'LarvaeSurvival_20μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt'] <- 'LarvaeSurvival_20μg-mLInsecticide_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'LarvaeSurvival_10μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt'] <- 'LarvaeSurvival_10μg-mLInsecticide_M'

'13-18'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt'] <- 'ElutionTime_Ethanol_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt'] <- 'ElutionTime_Ethanol_F'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt'] <- 'mCT_29C_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt'] <- 'mCT_26C_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt'] <- 'mCT_23C_M'
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'SolubleProteinLevels_PooledDiets_male.nogrms.original.txt'] <- 'ProteinLevels_PooledDiets_M'


ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt'] <- ''
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'Weight_PooledDiets_male.nogrms.original.txt'] <- ''
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'Weight_LowGlucoseDiet_male.nogrms.original.txt'] <- ''
ResampledResultsCulledExperimental$Phenotype[ResampledResultsCulledExperimental$Phenotype == 'Weight_HighGlucoseDiet_male.nogrms.original.txt'] <- ''


#boxplot
Plot_LeftBreakpoint <- ggplot(ResampledResultsCulledExperimental, aes(x=nTop_L_Breakpoint, y=Phenotype, color=factor(Iteration==1))) + 
  ggtitle("Left Breakpoint nTop") +
  labs(y = "Phenotype", x = "nTop") +
  geom_boxplot(alpha = 0.5)

#Plot_LeftBreakpoint <- ggplot(ResampledResultsCulledExperimental, aes(x=nTop_L_Breakpoint, y=Phenotype, color=factor(Iteration==1))) + 
  ggtitle("Left Breakpoint nTop") +
  geom_errorbar(aes(xmin=mean_nTop_LeftBreakpoint-sd_nTop_LeftBreakpoint, xmax = mean_nTop_LeftBreakpoint+sd_nTop_LeftBreakpoint), width=0.1, color="black", position="identity") +
  geom_point(alpha = 0.5) +
  labs(y = "Phenotype", x = "nTop")

#geom_errorbar(aes(ymin=mean_TotalFraction-sd, ymax=mean_TotalFraction+sd), width=0.1, color="black", position="identity") +
  


Plot_RightBreakpoint <- ggplot(ResampledResultsCulledExperimental, aes(x=nTop_R_Breakpoint, y=Phenotype, color=factor(Iteration==1))) + 
  ggtitle("Right Breakpoint nTop") +
  geom_boxplot(alpha = 0.5) +
  labs(y = "Phenotype", x = "nTop") 

Plot_INVMiddle <- ggplot(ResampledResultsCulledExperimental, aes(x=nTop_L_Breakpoint, y=Phenotype, color=factor(Iteration==1))) + 
  ggtitle("Inv Middle nTop") +
  geom_boxplot(alpha = 0.5) +
  labs(y = "Phenotype", x = "nTop") 
  

ggarrange(Plot_LeftBreakpoint, Plot_RightBreakpoint, Plot_INVMiddle + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)







'~~~~~~~~~~line graphs for enrichment aas a funciton of thr~~~~~~~~~~~~~~~~~~~~~'

VectorOf_2Lt_Pheno <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "ChillComaRecoveryTime_standard_female.nogrms.original.txt",
                        "ClimbingHeight_ParaquatInsecticideExposure_M.nogrms.jan.txt", "ClimbingHeight_standard_M.nogrms.jan.txt",
                        "DaySleep_standard_Female.nogrms.original.txt", 
                        "FoodIntake_standard_female.nogrms.original.txt", "FreeGlycerolLevels_PooledDiets_male.nogrms.original.txt", 
                        "FreeGlycerolLevels_LowGlucoseDiet_male.nogrms.original.txt", "FreeGlycerolLevels_HighGlucoseDiet_male.nogrms.original.txt", 
                        "LarvaeSurvival_40μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt", "LarvaeSurvival_20μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt", "LarvaeSurvival_10μg-mLChlorantraniliproleOn10mLFoodMedia_LarvaeMixed.nogrms.jan.txt",
                        "MeanElutionTime_RepeatedEthanolExposures_male.nogrms.original.txt", "MeanElutionTime_RepeatedEthanolExposures_female.nogrms.original.txt", 
                        "MinimumSurvivableTemperature_TwentyNineDegreesC_M.nogrms.jan.txt", "MinimumSurvivableTemperature_TwentySixDegreesC_M.nogrms.jan.txt", 
                        "MinimumSurvivableTemperature_TwentyThreeDegreesC_M.nogrms.jan.txt", "SolubleProteinLevels_PooledDiets_male.nogrms.original.txt", 
                        "SolubleProteinLevels_LowGlucoseDiet_male.nogrms.original.txt", "Weight_PooledDiets_male.nogrms.original.txt", 
                        "Weight_LowGlucoseDiet_male.nogrms.original.txt", "Weight_HighGlucoseDiet_male.nogrms.original.txt") 

VectorOf_2Lt_Pheno_Test <- c("ChillComaRecoveryTime_standard_male.nogrms.original.txt", "ChillComaRecoveryTime_standard_female.nogrms.original.txt")
ResultsEnrichmentSlidingThr <- data.table()

library(data.table)

# Create an empty data.table to store the results

for (j in seq_along(VectorOf_2Lt_Pheno)) {
  # Loop 100 times
  for (i in 1:6) {
    # Set the working directory and read in the data
    #i = 1
    #j = 1
    setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
    data <- fread(VectorOf_2Lt_Pheno[j])
    
    # Perform the calculations and save the results to a variable
    data[,rnp:=rank(PVAL)/(length(PVAL)+1)]
    data[,pos_bin:=round(POS, -4)]
    
    #left_inversion_boundary <- sample (1:16968732, 1)
    t <- if (i < 1.1 & i > 0.9) 0.05 else if (i < 2.1 & i > 1.9) 0.04 else if (i < 3.1 & i > 2.9) 0.03 else if (i < 4.1 & i > 3.9) 0.02 else if (i < 5.1 & i > 4.9) 0.01 else if (i < 6.1 & i > 5.9) 0.001
    
    #data[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
    data[CHR =="2L" & POS <= 2725744 & POS >= 1725744, inv:="Left Breakpoint"] #can also do these in thirds of the choromosome and it's also really itnersting
    data[CHR =="2L" & POS <= 12654180 & POS >= 2725744, inv:="2Lt-middle"]
    data[CHR =="2L" & POS <= 13654180 & POS >= 12654180, inv:="Right Breakpoint"]
    data[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
    data[, Threshold := t]
    #set(data, j="Threshold", value=Threshold)
    
    data2 <- data[, .(nTop = mean(rnp < Threshold), n = .N), keyby = .(CHR, inv)]
    data2[, thr := t]
    data2[, Iteration := i]
    data2[, p := pbinom(nTop * n, n, thr, lower.tail = F)]
    data2[, Phenotype := VectorOf_2Lt_Pheno[j]]
    data2[,thr:=t]
    data2[,PercentDifferenceFromThr:= ((nTop - thr)/((nTop + thr)/2))*100] #can add abs() for top of fraction for positive percent difference
    
    # Append the results to the 'ResultsEnrichmentSlidingThr' data.table
    ResultsEnrichmentSlidingThr <- rbindlist(list(ResultsEnrichmentSlidingThr, data2))
  }
}
write.csv(ResultsEnrichmentSlidingThr,"/Users/dgg/backupForEnrichmentAcrossThr.csv")
ResultsEnrichmentSlidingThr <- fread("backupForEnrichmentAcrossThr.csv")
ResultsEnrichmentSlidingThr <- ResultsEnrichmentSlidingThr[complete.cases(ResultsEnrichmentSlidingThr$inv)]



ggplot(ResultsEnrichmentSlidingThr, aes(x=thr, y=nTop, color=inv)) + 
  ggtitle("Contained Fraction of Total Significant SNPs over Varying Significance Thresholds") +
  geom_point(alpha=0.5) +
  geom_line(aes(group=inv)) +  # add a line for each 'inv' type
  scale_color_manual(values=c("blue","purple","red","green"), labels=c("2Lt-left edge","2Lt-middle","2Lt-right edge","outside-2Lt")) + # specify colors and labels for the key
  theme(legend.position="top") +
  facet_wrap(~Phenotype) 


ggplot(ResultsEnrichmentSlidingThr, aes(x=thr, y=PercentDifferenceFromThr, color=inv)) + 
  facet_wrap(~Phenotype) +
  ggtitle("Percent Difference between contained Total Significant SNPs to Null Expectation") +
  labs(y = "Percent Difference Between Observed Top SNP Concentration and Null Expectation", x = "Null Expectation Thresholds (0.001, 0.01, 0.02, 0.03, 0.04, 0.05)") +
  geom_point(alpha=0.5) +
  geom_line(aes(group=inv)) +  # add a line for each 'inv' type
  scale_color_manual(values=c("blue","purple","red","green"), labels=c("Left Breakpoint","2Lt-middle","Right Breakpoint","outside-2Lt")) + # specify colors and labels for the key
  theme(legend.position="top") +
  ylim(-125, 125) + 
  geom_hline(yintercept = 0)


ResultsEnrichmentSlidingThrGrouped <- group_by(ResultsEnrichmentSlidingThr, inv)
ResultsEnrichmentSlidingThrGroupedsummary <- summarize(ResultsEnrichmentSlidingThrGrouped, mean_nTop = mean(nTop))


ResultsEnrichmentSlidingThr %>%
  group_by(inv) %>%
  summarise(avg = mean(nTop))









