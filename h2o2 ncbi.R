# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Dec 3 22:51:59 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
# Read MicroArray data from local 

gset <- getGEO(file="GSE19983_series_matrix.txt.gz", GSEMatrix =TRUE, AnnotGPL=FALSE)

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for Control and H2O2 treat 1 hr
gsms <- "XXX111XXX000XXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", number=42417)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SEQUENCE","GB_ACC","ORF","miRNA_ID"))

# Clean junk variable
rm(fl,sml,sel,LogC,qx,gsms,i,cont.matrix,design,fit,fit2)

##################################################################################################
library(magrittr)

# select differential express gene in MicroArray
# True : differential gene
select <- tT$adj.P.Val<0.05 & (tT$logFC>1|tT$logFC< -1)

# -log10(adj.pvalue) transformation
LOG10 <- sapply(tT$adj.P.Val,function(x){
  return(-log10(x))
})

tT <- cbind(tT,LOG10)
tT <- cbind(tT,select)

# Fillter by select
gene <- tT[tT$select,]

# Extract RAP-ID from gene$ID
name <- sapply(gene$ID,function(x){
  ID = strsplit(x,split="|",fixed=T)
  ID = ID[[1]][1]
})
name <- data.frame(name)
gene <- cbind(name,gene)

########################################################################################

## Visualization of MicroArray data
library(ggplot2)
qplot(tT$logFC,tT$adj.P.Val)
ggplot(tT,aes(x=logFC,y=LOG10,color=select),xlim=c(-2.5,2)) + geom_point()
t$logFC
write.table(tT, file="H2O2tT.csv", row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO(file="GSE19983_series_matrix.txt.gz", GSEMatrix =TRUE, AnnotGPL=FALSE)


# group names for all samples in a series
gsms <- "XXX111XXX000XXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("cont","h2o2")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE19983", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")


