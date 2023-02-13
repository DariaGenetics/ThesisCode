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
  data <- fread(VectorOfPheno[i])


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

  data2[,Phentope:= VectorOfPheno[i]] #can use a split to take off file extensions from names
  return(data2)
}
gather = rbindlist(PhenoGather)


ggplot(data=gather, aes(x=pos_bin, y=-log10(p), color=inv)) + geom_line() + facet_grid(Phentope~CHR)
