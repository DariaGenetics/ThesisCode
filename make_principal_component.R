### libraries
  library(missMDA)
  library(ggplot2)
  library(data.table)
  library(FactoMineR)
  library(factoextra)

### load data
  setwd("/Users/alanbergland/Documents/GitHub/ThesisCode")
  AllPhenoWInversionStatus <- fread("AllPheno.csv", fill = TRUE)
  AllPhenoWInversionStatus <- AllPhenoWInversionStatus[c(1:201),]

### clean up
  AllPhenoWInversionStatus <- subset(AllPhenoWInversionStatus, select = -c(V233, V234, V235))
  RalID <- as.vector(unlist(AllPhenoWInversionStatus$ral_id))
  In2LtStatus <- as.vector(unlist(AllPhenoWInversionStatus$In2Lt))

  AllPhenoWOInversionStatus = subset(AllPhenoWInversionStatus, select = -c(In2Lt, ral_id) ) #colnumber [231]
  AllPhenoWOInversionStatus <- as.matrix(AllPhenoWOInversionStatus)
  dimnames(AllPhenoWOInversionStatus)[[1]] <- paste("RAL_", RalID, sep="")

### impute
  AllPhenoPCA1 = imputePCA(AllPhenoWOInversionStatus, scale.unit=TRUE, ncp=5, graph=T)

### check - is the output of this function the PC or the imputed data?
  str(AllPhenoPCA1)
  table(apply(AllPhenoWOInversionStatus, 1, function(x) mean(is.na(x))))


  apply(AllPhenoPCA1$completeObs, 2, function(x) mean(is.na(x)))
  table(apply(AllPhenoPCA1$fittedX, 2, function(x) mean(is.na(x))))

  plot(AllPhenoPCA1$completeObs[,200] ~ AllPhenoPCA1$fittedX[,200])
  plot(AllPhenoWOInversionStatus[,200]$StarvationResistance_standard_male ~ AllPhenoPCA1$completeObs[,200])
  plot(AllPhenoWOInversionStatus[,200]$StarvationResistance_standard_male ~ AllPhenoPCA1$fittedX[,200])

###run PCA
  pc <- PCA(AllPhenoPCA1, ncp=250)
  dim(pc$ind$coord)
  str(pc$ind$coord)

### collect output and merge with inversion status
  pc_out <- as.data.table(pc$ind$coord)
  pc_out[,ral_id:=row.names(pc$ind$coord)]
  inv <- AllPhenoWInversionStatus[,c("ral_id", "In2Lt"), with=F]
  inv[,ral_id:=paste("RAL_", ral_id, sep="")]

  pc_out <- merge(pc_out, inv, by="ral_id")

  plot.PCA(pc)

### save
  save(pc_out, pc, file="/Users/alanbergland/Documents/GitHub/ThesisCode/pca_output.Rdata")


  PrincipleComponents <- AllPhenoPCA1$fittedX
PrincipleComponents2 <- as.data.frame(PrincipleComponents) #note - V1 = principle component 1, V2 = principle component 2, etc...
PrincipleComponents2$RalId <- RalID
PrincipleComponents2$In2LtStatus <- In2LtStatus
