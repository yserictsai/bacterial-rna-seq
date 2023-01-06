#!/bin/bash



-----------
#FastQC check
 

for i in {1..6};do
    fastqc -o ../fastqc MYQT${i}_Clean_Data1.fq.gz
    fastqc -o ../fastqc MYQT${i}_Clean_Data2.fq.gz
done

#Note: 


----------

#Trimmomatics preprocessing
# trim adaptor (這個dataset 不需要)
# filter low complexity reads (這個dataset不需要)
# trim polyA-tails (原核沒有polyA tails)
# GC content (can't remove by Preprocessing, 原核會稍微偏高)
# [remove E coli. sequences, sequences contamination] (不確定需不需要) 



#Trimmomatics
java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog QC_Preprocessing.log MYQT1_Clean_Data1.fq MYQT1_Clean_Data2.fq paired_MYQT1_Clean_Data1.fq unpaired_MYQT1_Clean_Data1.fq paired_MYQT1_Clean_Data2.fq unpaired_MYQT1_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25


#Note:
# remove low quality base reads (去除Q-score平均低於25的reads), read average quality < 25 filter (AVGQUAL:25)
# Read length control (短於50bp的read去除), filter length < 50bp read (MINLEN:50)
# sequence specific bias caused by random hexamer (去除5'end的第一個base以校正), remove the first 5' base (HEADCROP:1)
# 3' end when the base quality < 25 trim (TRAILING:25)
# 3' start when the base quality < 25 trim (LEADING:25)
# slidingwindow 3-based < 20 mean read quality (SLIDINGWINDOW:3:20)


#PRINSEQ

perl /app/programs/prinseq-lite-0.20.4/prinseq-lite.pl -fastq paired_MYQT1_Clean_Data1.fq -fastq2 paired_MYQT1_Clean_Data2.fq -ns_max_n 2 -derep 1 -derep_min 150  -out_good PRINSEQfiltered -out_bad null -no_qual_header -log -verbose

# remove PCR artifacts (去除發生150次以上的duplicates 不去除太低的以避免破壞gene expression的dynamics), (-derep 1 -derep_min 150)
# filter the read has Ambiguous Base (去除同一條read出現兩個以上Ns base), (-ns_max_n 2)


----------

#Bowtie2 mapping

bowtie2-build ./reference/SM_db11.fa ./reference/SM_db11

bowtie2 -p 8 -x ./reference/SM_db11 -1 ./fastq/control/paired_MYQT1_Clean_Data1.fq -2 ./fastq/control/paired_MYQT1_Clean_Data2.fq -S ./bam/MYQT1.bowtie2.sam && samtools view -bS ./bam/MYQT1.bowtie2.sam > ./bam/MYQT1.bowtie2.bam  && rm ./bam/MYQT1.bowtie2.sam

#Note:


----------
#view through IGV


#Note



----------
#eXpress Quantification

for i in {1..6};do 
	express --output-dir ./express/MYQT${i}_express ./reference/SM_db11.fa ./bam/MYQT${i}.bowtie2.bam &  
done


#Note: