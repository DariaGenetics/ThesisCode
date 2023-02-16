setwd("/Users/alanbergland/Daria Pheno")

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)

#CentroidSize_Standard_F <- read.delim("CentroidSize_standard_F.nogrms.jan.txt",
#                                      stringsAsFactor = FALSE)
#WingCentroidSize_standard_F <- read.delim("WingCentroidSize_standard_female.nogrms.original.txt",
#                                          stringsAsFactor = FALSE)
#
#attempting a gather

#VectorOfPheno<- c("CentroidSize_standard_F.nogrms.jan.txt", "WingCentroidSize_standard_female.nogrms.original.txt") #add the rest of the file names here
VectorOfPheno<- list.files("/Users/alanbergland/Daria Pheno")

PhenoGather <- foreach(i=1:length(VectorOfPheno)) %dopar% {

  #i = 1
  #setwd("/Users/dgg/Desktop/Thesis/GWAS STUFF/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes")
  setwd("/Users/alanbergland/Daria Pheno")


  #data = VectorOfPheno[i] #set data
  #data = as.data.table(read.delim(data)) #loads data
  data1 <- fread(VectorOfPheno[1])
  data1[,rnp_1:=rank(PVAL)/(length(PVAL)+1)]
  data1[,pos_bin:=round(POS, -4)]

  data1[CHR =="2L" & POS <= 13154180 & POS >= 2225744, inv:="2Lt"]
  data1[CHR =="2L" & is.na(inv), inv:="outside-2Lt"]
  data1[CHR =="3R" & POS <= 29031297 & POS >= 21406917, inv:="3RMo"]
  data1[CHR =="3R" & is.na(inv), inv:="outside-3RMo"]

  data2 <- fread(VectorOfPheno[2])
  data2[,rnp_2:=rank(PVAL)/(length(PVAL)+1)]
  data2[,pos_bin:=round(POS, -4)]

  setkey(data1, "CHR", "POS")
  setkey(data2, "CHR", "POS")

  data <- merge(data1, data2)

  data2 <- data[,list(nTop = mean(rnp_1 < 0.05 & rnp_2 < .05), thr = c(0.05^2), n=length(rnp_1), eq=mean(rnp_1==rnp_2)),
                list(CHR, pos_bin=pos_bin.x, inv)]

  data2[,p:=pbinom(nTop*n, n, thr, lower.tail=F)]

  data2[,Phentope:= VectorOfPheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
gather = rbindlist(PhenoGather)


w <- dcast()

ggplot(data=data2, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(.~CHR)
