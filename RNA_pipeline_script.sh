#!/usr/bin/bash

# test fastq-dump for 5 lines
fastq-dump -X 5 -Z --skip-technical  --readids --dumpbase  --clip SRR1746771


# fastq-dump the dataset

fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746771
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746772
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746773
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746774
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746775
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746776
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746777
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746778
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746779
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746780
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746781
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746782
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746783
fastq-dump --outdir ./fastq --gzip --skip-technical  --readids --dumpbase  --clip SRR1746784


#fastQC
chmod 755 fastqc

-a adaptor_file [name(tab)sequence]
-c contaminants_file [name(tab)sequence]

fastqc -o ../fastqc SRR1767248.fastq.gz
fastqc -o ../fastqc SRR1767249.fastq.gz
fastqc -o ../fastqc SRR1767250.fastq.gz
fastqc -o ../fastqc SRR1767251.fastq.gz
fastqc -o ../fastqc SRR1767252.fastq.gz
fastqc -o ../fastqc SRR1767253.fastq.gz
fastqc -o ../fastqc SRR1767254.fastq.gz
fastqc -o ../fastqc SRR1767255.fastq.gz
fastqc -o ../fastqc SRR1767256.fastq.gz
fastqc -o ../fastqc SRR1767257.fastq.gz
fastqc -o ../fastqc SRR1767258.fastq.gz
fastqc -o ../fastqc SRR1767259.fastq.gz
fastqc -o ../fastqc SRR1767260.fastq.gz
fastqc -o ../fastqc SRR1767261.fastq.gz
fastqc -o ../fastqc SRR1767262.fastq.gz
fastqc -o ../fastqc SRR1767263.fastq.gz

#prinseq

要裝JSON.pm, Cairo.pm

perl ~/tools/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ../fastq/SRR1767244.fastq -out_good null -out_bad null -graph_data graph
perl ~/tools/prinseq-lite-0.20.4/prinseq-graphs.pl -i graph -html_all -o QCreport

#trimmomatic - total
#trim - filter by length - filter read mean quality - sequence specific bias - Ambiguous base - remove duplicates - sequence contamination - low complexity
java -jar /home/ystsai/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 5 -phred33 -trimlog QC.log ../fastq/SRR1767244.fastq.gz QC_SRR1767244.fastq.gz HEADCROP:1 SLIDINGWINDOW:3:22 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25

#PRINSEQ - total 不能吃gz

perl ~/tools/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ../trimmomatic/QC_SRR1767244.fastq -ns_max_n 2 -derep 1 -derep_min 150 -lc_threshold 7 -lc_method dust -out_good PRINSEQfilter_SRR1767244 -out_bad null -no_qual_header -log -verbose

#low complexity
perl ~/tools/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ../trimmomatic/QC_SRR1767244.fastq -ns_max_n 2 -derep 1 -derep_min 101 -lc_threshold 7 -lc_method dust -out_good PRINSEQfilter_SRR1767244 -out_bad null -no_qual_header -log -verbose

#PolyA

#FastQ



##Base quality filter for  SE - read average quality < 25 filter
java -jar /home/ystsai/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 -trimlog QC.log SRR1767244.fastq.gz OUT_SRR1767244.fastq.gz AVGQUAL:25

##Trimming for SE - 3' end when the base quality < 25 trim (TRAILING:25), 3' start when the base quality < 25 trim (LEADING:25), filter length < 50bp read (MINLEN:50), slidingwindow 3-based < 20 mean read quality .
java -jar /home/ystsai/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 -trimlog QC.log SRR1767244.fastq.gz OUT_SRR1767244.fastq.gz


##Ambiguous base - remove the read has more than 2 Ns  fastQC 會檢查 或沒過 再執行也可
prinseq-lite.pl -fastq ../trimmomatic/Trim_SRR1767244.fastq.gz -ns_max_n 2 -out_good nfilter_SRR1767244 -out_bad null -no_qual_header -log -verbose

## remove adaptor 如果fastqc pass 就不用做 若沒過fastqc 會說要用哪個fa
java -jar /home/ystsai/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 -trimlog QC.log SRR1767244.fastq.gz OUT_SRR1767244.fastq.gz ILLUMINACLIP:/home/ystsai/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10

## bias caused by random hexamer priming - 用trimmomatic 去除第一個base

java -jar /home/ystsai/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 -trimlog QC.log SRR1767244.fastq.gz OUT_SRR1767244.fastq.gz HEADCROP:1

## remove duplicates - removes exact duplicate reads > 100 times
prinseq-lite.pl -fastq ../fastq/SRR1767244.fastq -derep 1 -derep_min 101 -log -verbose -out_good dupfilter_SRR1767244 -out_bad null -no_qual_header

## remove low-complexity  DUST >7   ENTROPY <50  polyA tail -trim_tail_right 5  依照PRINSEQ 的報告決定

prinseq-lite.pl -fastq ../fastq/SRR1767244.fastq -lc-threshold 50 -lc-method entropy -log -verbose -out_good dupfilter_SRR1767244 -out_bad null -no_qual_header
or
prinseq-lite.pl -fastq ../fastq/SRR1767244.fastq -lc-threshold 7 -lc-method dust -log -verbose -out_good dupfilter_SRR1767244 -out_bad null -no_qual_header


###############
trim - filter by length - filter read mean quality - sequence specific bias -  remove duplicates - sequence contamination - low complexity
