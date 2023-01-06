#https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#https://cgrlucb.wikispaces.com/edgeRWikiVersion
#https://pods.iplantcollaborative.org/wiki/pages/viewpage.action?pageId=6720810
#https://github.com/nibb-gitc/gitc2015a-rnaseq/wiki/ex5

#DESeq vs edgeR Comparison
#http://www.r-bloggers.com/deseq-vs-edger-comparison/

library("edgeR")
library("DESeq2")
library("ggplot2")
setwd("~/Google Drive/bugInHost.Mazmanian/results/032216_V3data/bug")
#get functional annotations
functionalAnnotation<-read.csv("/Users/wcchou/Google Drive/bugInHost.Mazmanian/data/reference/finalLatOriIntAnnotation.csv", header=T)

# work on v3 RNA-seq data
directory <- "./"

# work on v1's wild-type RNA-seq data

# Lumen
(sampleFiles <- grep("counts", list.files(directory),value=T))
(sampleFiles <- grep("HS",sampleFiles, value=T))
(sampleFiles <- grep("lumen",sampleFiles, value=T))
# Mucus
(sampleFiles <- grep("counts", list.files(directory),value=T))
(sampleFiles <- grep("HS",sampleFiles, value=T))
(sampleFiles <- grep("mucus",sampleFiles, value=T))
# Tissue
(sampleFiles <- grep("counts", list.files(directory),value=T))
(sampleFiles <- grep("HS",sampleFiles, value=T))
(sampleFiles <- grep("tissue",sampleFiles, value=T))

listALL=NULL
for(i in 1:length(sampleFiles)){
    x=read.table(sampleFiles[i],header=F,sep="\t")
    listALL[[i]]=x[,2]
}
length(listALL)

WTvsCCF <- do.call(cbind, listALL)
rownames(WTvsCCF)=x[,1]
x=WTvsCCF
#group <- factor(c(1,1,1,2,2,2))
group <- c(rep("WT",3),rep("CCF",3))
y <- DGEList(counts=x,group=group)
design <- model.matrix(~group)
dim(y)
#keep only those genes that have at least 1 read per million in at least 3 samples.
# y <- y[rowSums(1e+06 * y$counts/expandAsMatrix(y$samples$lib.size, dim(y)) > 1) >= 3, ]
# dim(y)
(y <- calcNormFactors(y))

(y <- estimateCommonDisp(y))
(y <- estimateTrendedDisp(y))
(y <- estimateTagwiseDisp(y))

# (y <- estimateCommonDisp(y, design))
# (y <- estimateTrendedDisp(y, design))
# (y <- estimateTagwiseDisp(y, design))

#(y <- estimateDisp(y))
#(y <- estimateDisp(y,design))

# y <- estimateGLMCommonDisp(y, design)
# y <- estimateGLMTrendedDisp(y, design) 
# y <- estimateGLMTagwiseDisp(y, design)
# 
# y <- estimateGLMCommonDisp(y)
# y <- estimateGLMTrendedDisp(y) 
# y <- estimateGLMTagwiseDisp(y)


# y <- estimateGLMTrendedDisp(y, design )
# y <- estimateTagwiseDisp(y, design )
#y <- estimateCommonDisp(y,design)

plotBCV(y, cex=0.5, main="edgeR: Biological coefficient of variation (BCV) vs abundance")


# #glmQLFit
# fit <- glmQLFit(y,design)
# qlf <- glmQLFTest(fit,coef=2)
# #topTags(qlf, n=20)
# eqtable <- topTags(qlf, n=nrow(qlf))$table
# eqtable <- eqtable[order(eqtable$FDR), ]
# head(eqtable)
# dim(eqtable)
# dim(subset(eqtable, FDR<0.05))
# subset(eqtable, FDR<0.05)
# etable <- eqtable

#glmFit
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt, n=20)
## Make a table of results
etable <- topTags(lrt, n=nrow(lrt))$table
normalizedCount<-y$pseudo.counts
colnames(normalizedCount)=sub(".bug.bowtie2.psorted.coverage.counts","",sampleFiles)
head(normalizedCount)
etable <- merge(etable,normalizedCount, by='row.names' )
etable <- etable[order(etable$FDR), ]
#dim(etable)
etable=etable[-which(duplicated(etable$Row.names)),]
row.names(etable)=etable$Row.names
etable=etable[,-1]
head(etable,20)
dim(etable)
dim(subset(etable, FDR<=0.05))

subset_etable=merge(subset(etable, FDR<=0.05), functionalAnnotation, by.x="row.names", by.y="latestLocusTag")
colnames(subset_etable)[1]="locusTag"
subset_etable <- subset_etable[order(subset_etable$FDR), ]

#write.csv(subset_etable,"./Bacteria_WTvsCCF_Lumen.edgeR.csv", row.names = F)
