### libraries
  library(data.table)
  library(effectsize)
  library(ggplot2)
  library(foreach)
  library(boot)
  library(doMC)
  registerDoMC(4)

### load data
  load("/Users/alanbergland/Documents/GitHub/ThesisCode/pca_output.Rdata")

  PrincipleComponentsTrimmed <- fread("/Users/alanbergland/Documents/GitHub/ThesisCode/PrincipleComponentsWOLastFourRows.csv")
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

### which factors contribute to Dim.8?
  loading <- as.data.table(pc$var$coord)
  loading[,pheno:=row.names(pc$var$coord)]
  ll <- melt(loading, id.var="pheno")

  ll.rank <- ll[,list(quan=rank(abs(value))/(length(value)+1), pheno=pheno, rank=rank(value)), list(variable)]

  setkey(ll, variable, pheno)
  setkey(ll.rank, variable, pheno)
  ll <- merge(ll, ll.rank)

  phenoLoad <- ggplot(data=ll[variable=="Dim.8"][quan>.9]) +
  geom_segment(aes(y=rank(rank), yend=rank(rank), x=0, xend=value)) +
  geom_text(data=ll[variable=="Dim.8"][quan>.9],
              aes(y=rank(rank), x=.5, label=pheno), hjust=0, size=2) +
  xlim(c(-.5, 2.5)) + xlab("Loading value") + ylab("phenotype")



### effect plot
  pc.ag <- MeltedPrincipleComponents[, list(mu=mean(ComponentValues), sd=sd(ComponentValues), n=length(ComponentValues)),
                              list(In2Lt, PC)]
  pc.ag[,se:=sd/sqrt(n)]

  ggplot(data=pc.ag[PC=="Dim.8"], aes(x=In2Lt, y=mu)) +
  geom_point() +
  geom_segment(aes(x=In2Lt, xend=In2Lt, y=mu-2*se, yend=mu+2*se))



ggplot(data=A, aes(x=PC, y=r2)) +
geom_segment(aes(x=PC, xend=PC, y=uci, yend=lci)) +
geom_point() +
geom_point(data=A[p<.05], aes(x=PC, y=r2), color="red") +
geom_point(data=A[pa<.1], aes(x=PC, y=r2), color="blue")








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
