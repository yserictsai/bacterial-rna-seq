# Strategy
#
# 1. Genome-mapping (gapped mapper)
# Reads - map to genome reference (Tophat, STAR) - with GFF then Transcript identification & counting.
# Reads - map to genome reference (Tophat, STAR) - without GFF then Transcript discovery & counting. - Blast2GO (functional annotation)
#
# 2. Genome-mapping (gapped mapper) and Genome-guided Transcriptome assembly
# Reads - map to genome reference (Tophat, STAR) - assembly (cufflink, Scripture) - with GFF then Transcript identification & counting.
# Reads - map to genome reference (Tophat, STAR) - assembly (cufflink, Scripture) - without GFF then Transcript discovery & counting. - Blast2GO (functional annotation)
#
#
# 3. Transcriptome-mapping （ungapped mapper)
# Reads - map to Transcriptome reference (cds, cdna)(Bowtie) -  Transcript identification & counting.
#
#
# 4. De novo assembly
# Reads - Assembly to Transcript - map reads back (bowtie) - GTF-based counting - Blast2GO (functional annotation)

# Note
# illumina 1.8+, sanger = 33; illumina 1.3+ 1.5+ = 64




#Step 1 QC-Check

#!/bin/bash

for i in {1..6};do
    fastqc -o ../fastqc MYQT${i}_Clean_Data1.fq.gz
    fastqc -o ../fastqc MYQT${i}_Clean_Data2.fq.gz
done


perl ~/tools/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ../fastq/SRR1767244.fastq -out_good null -out_bad null -graph_data graph
perl ~/tools/prinseq-lite-0.20.4/prinseq-graphs.pl -i graph -html_all -o QCreport



#Step 2 Data Preprocessing


filtering - trimming - error correction - bias correction

# remove low quality base reads (去除Q-score平均低於25的reads)
# filter the read has Ambiguous Base (去除同一條read出現兩個以上Ns base)
# trim adaptor (這個dataset 不需要)
# Read length control (短於50bp的read去除)
# sequence specific bias caused by random hexamer (去除5'end的第一個base以校正)
# remove PCR artifacts (去除發生150次以上的duplicates 不去除太低的以避免破壞gene expression的dynamics)
# [remove E coli. sequences, sequences contamination] (不確定需不需要)
# filter low complexity reads (不確定需不需要)
# trim polyA-tails (不確定需不需要)
# GC content (can't remove by Preprocessing)
# read average quality < 25 filter (AVGQUAL:25)
# 3' end when the base quality < 25 trim (TRAILING:25)
# 3' start when the base quality < 25 trim (LEADING:25)
# filter length < 50bp read (MINLEN:50)
# slidingwindow 3-based < 20 mean read quality (SLIDINGWINDOW:3:20)
# remove the first 5' base (HEADCROP:1)
	java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog MYQT4_QC_Preprocessing.log ./fastq/treatment/MYQT4_Clean_Data1.fq.gz ./fastq/treatment/MYQT4_Clean_Data2.fq.gz paired_MYQT4_Clean_Data1.fq unpaired_MYQT4_Clean_Data1.fq paired_MYQT4_Clean_Data2.fq unpaired_MYQT4_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25


for i in {4..6};do

	java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog MYQT${i}_QC_Preprocessing.log ./fastq/treatment/MYQT${i}_Clean_Data1.fq.gz ./fastq/treatment/MYQT${i}_Clean_Data2.fq.gz paired_MYQT${i}_Clean_Data1.fq unpaired_MYQT${i}_Clean_Data1.fq paired_MYQT${i}_Clean_Data2.fq unpaired_MYQT${i}_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25


done

for i in {4..6};do

	java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog MYQT${i}_QC_Preprocessing.log MYQT${i}_Clean_Data1.fq MYQT${i}_Clean_Data2.fq paired_MYQT${i}_Clean_Data1.fq unpaired_MYQT${i}_Clean_Data1.fq paired_MYQT${i}_Clean_Data2.fq unpaired_MYQT${i}_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25


done

# filter out the reads which have more than two Ns (-ns_max_n 2)
# if the reads have adaptor (ILLUMINACLIP)
# remove exact duplicate reads which occur more than 150 times (-derep 1 -derep_min 151)
# filter low complexity reads (-lc_threshold 7 -lc_method dust or -lc_threshold 50 -lc_method entropy)
# trim polyA-tails (-trim_tail_right or -trim_tail_left)

for i in {1..3};do

	perl ~/tools/prinseq-lite-0.20.4/prinseq-lite.pl -fastq paired_MYQT${i}_Clean_Data1.fq -fastq2 paired_MYQT${i}_Clean_Data2.fq \
	-ns_max_n 2 -derep 1 -derep_min 150 -lc_threshold 7 -lc_method dust -out_good PRINSEQfiltered -out_bad null -no_qual_header -log -verbose

done



	for i in {1..6};do
		java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog QC_Preprocessing.log MYQT${i}_Clean_Data1.fq MYQT${i}_Clean_Data2.fq paired_MYQT${i}_Clean_Data1.fq unpaired_MYQT${i}_Clean_Data1.fq paired_MYQT${i}_Clean_Data2.fq unpaired_MYQT${i}_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25
	done


Mapping using bowtie2

#bowtie2-build [options]* <reference_in> <bt2_base>




java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog QC_Preprocessing.log MYQT1_Clean_Data1.fq MYQT1_Clean_Data2.fq paired_MYQT1_Clean_Data1.fq unpaired_MYQT1_Clean_Data1.fq paired_MYQT1_Clean_Data2.fq unpaired_MYQT1_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25


perl /app/programs/prinseq-lite-0.20.4/prinseq-lite.pl -fastq paired_MYQT1_Clean_Data1.fq -fastq2 paired_MYQT1_Clean_Data2.fq -ns_max_n 2 -derep 1 -derep_min 150  -out_good PRINSEQfiltered -out_bad null -no_qual_header -log -verbose



bowtie2-build ./reference/SM_db11.fa ./reference/SM_db11

bowtie2 -p 8 -x ./reference/SM_db11 -1 ./fastq/control/paired_MYQT1_Clean_Data1.fq -2 ./fastq/control/paired_MYQT1_Clean_Data2.fq -S ./bam/MYQT1.bowtie2.sam && samtools view -bS ./bam/MYQT1.bowtie2.sam > ./bam/MYQT1.bowtie2.bam  && rm ./bam/MYQT1.bowtie2.sam


for i in {1..3};do

bowtie2 -p 8 -x ./reference/SM_db11 -1 ./fastq/control/MYQT${i}_Clean_Data1.fq -2 ./fastq/control/MYQT${i}_Clean_Data2.fq -S ./bam/MYQT${i}.bowtie2.sam && samtools view -bS ./bam/MYQT${i}.bowtie2.sam > ./bam/MYQT${i}.bowtie2.bam  && rm ./bam/MYQT${i}.bowtie2.sam

done



for i in {4..6};do

bowtie2 -p 8 -x ./reference/SM_db11 -1 ./fastq/treatment/MYQT${i}_Clean_Data1.fq -2 ./fastq/treatment/MYQT${i}_Clean_Data2.fq -S ./bam/MYQT${i}.bowtie2.sam && samtools view -bS ./bam/MYQT${i}.bowtie2.sam > ./bam/MYQT${i}.bowtie2.bam  && rm ./bam/MYQT${i}.bowtie2.sam

done





#Post-mapping QC


RseQC



#Quantification using bedtools


Gene_id  start  end 











# DEG using edgeR

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




Post mapping quality control

coverage uniformity along transcripts

saturation of sequencing depth

ribosomal RNA content

reads counted per genes, sample relations, batch effect -> Heatmaps, PCA plots.
