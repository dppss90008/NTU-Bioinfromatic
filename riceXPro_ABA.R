# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Dec 6 02:14:01 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(magrittr)

# load series and platform data from GEO

gset <- getGEO(file="GSE39429_series_matrix.txt.gz", GSEMatrix =TRUE, AnnotGPL=FALSE)
#if (length(gset) > 1) idx <- grep("GPL6864", attr(gset, "names")) else idx <- 1
#gset <- gset[[idx]]
gset2 <- gset
gset <- gset2

Hormone <- gset$description %>% as.vector()
gsms <- c()
for (i in Hormone) {
  if (grepl("root at 0 min",i)) {
    gsms <- c(gsms,"0")
  }else if (grepl("root at 1 hr after jasmonic acid",i)){
    gsms <- c(gsms,"1")
  }else{
    gsms <- c(gsms,"X")
  }
}

sml <- gsms
#gsms

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
#gsms <- paste0("000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
#               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
#               "XXXXXXXXXXXXXXXXXXXXXXXXXX111XXXXXXXXXXXXXX")
#sml <- c()
#for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

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
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(ex))

tT <- tT[tT$SEQUENCE!="",]


select <- tT$P.Value<0.05 & (tT$logFC>1|tT$logFC< -1)
LOG10 <- sapply(tT$adj.P.Val,function(x){
  return(-log10(x))
})

tT <- cbind(tT,LOG10)
tT <- cbind(tT,select)
gene <- tT[tT$select,]

library(magrittr)
name <- data.frame(strsplit(gene$Accessions,split="|",fixed=T)) %>% t()
gene <- cbind(gene,name[,1])

name <- sapply(gene$Accessions,function(x){
  ID = strsplit(x,split="|",fixed=T)
  ID = ID[[1]][1]
})

name <- data.frame(name)
write.csv(name,file="JA_Shoot.csv")
  
library(ggplot2)

ggplot(tT,aes(x=logFC,y=LOG10,color=select),xlim=c(-2.5,2)) + geom_point() + ylab("-log10(P.value)")




