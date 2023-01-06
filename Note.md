Bacteria RNA-seq Pipeline


github: 
	
	https://github.com/wenchichou/bacterialRNAseq.git

ftp:
	
	ftp.broadinstitute.org

Raw data:
	
	serratia marcescens db11, 3 x Control (MYQT{1..3}), 3 x Treatment (MYQT{4..6}).

Reference:
	
	Source link:

	GCF_000513215.1_DB11_genomic.fna

Annotation: 
	
	Source link:

	GCF_000513215.1_DB11_genomic.gff

----------------------------------------------------------------------------
Step1 FastQC check

Input:
	
	MYQT${i}_Clean_Data1.fq.gz
	MYQT${i}_Clean_Data2.fq.gz

Output:

	MYQT${i}_Clean_Data1_fastqc.html

SourceCode:

	for i in {1..6};do
    	fastqc -o ../fastqc MYQT${i}_Clean_Data1.fq.gz
    	fastqc -o ../fastqc MYQT${i}_Clean_Data2.fq.gz
	done


Result:

	preFastQC_result.zip


Note: 



----------------------------------------------------------------------------
Step2 Trimmomatics preprocessing


Input:
	
	MYQT${i}_Clean_Data1.fq.gz
	MYQT${i}_Clean_Data2.fq.gz

Output:

	QC_MYQT${i}_Preprocessing.log 
	paired_MYQT${i}_Clean_Data1.fq 
	paired_MYQT${i}_Clean_Data2.fq 
	unpaired_MYQT${i}_Clean_Data1.fq 
	unpaired_MYQT${i}_Clean_Data2.fq

SourceCode:

	for i in {1..6};do
		java -jar /app/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog QC_Preprocessing.log MYQT${i}_Clean_Data1.fq MYQT${i}_Clean_Data2.fq paired_MYQT${i}_Clean_Data1.fq unpaired_MYQT${i}_Clean_Data1.fq paired_MYQT${i}_Clean_Data2.fq unpaired_MYQT${i}_Clean_Data2.fq HEADCROP:1 SLIDINGWINDOW:3:20 LEADING:22 TRAILING:22 MINLEN:50 AVGQUAL:25
	done

Result:

	xxx.zip  


Note: 

	trim adaptor (原始資料已完成) 
	filter low complexity reads (不需要)
	trim polyA-tails (原核沒有polyA tails)     
	GC content (Can't remove by Preprocessing, 原核會稍微偏高)
	remove sequences contamination (不確定需不需要) 
	remove low quality base reads (去除Q-score平均低於25的reads), read average quality < 25 filter (AVGQUAL:25)
	Read length control (短於50bp的read去除), filter length < 50bp read (MINLEN:50)
	sequence specific bias caused by random hexamer (去除5'end的第一個base以校正), remove the first 5' base (HEADCROP:1)
	3' end when the base quality < 25 trim (TRAILING:25)
	3' start when the base quality < 25 trim (LEADING:25)
	slidingwindow 3-based < 20 mean read quality (SLIDINGWINDOW:3:20)
	remove PCR artifacts (去除發生150次以上的duplicates 不去除太低的以避免破壞gene expression的dynamics), (-derep 1 -derep_min 150) (不確定需不需要)
	filter the read has Ambiguous Base (去除同一條read出現兩個以上Ns base), (-ns_max_n 2) (不確定需不需要)


----------------------------------------------------------------------------

Step3 Post-FastQC


Input:
	
	paired_MYQT${i}_Clean_Data1.fq 
	paired_MYQT${i}_Clean_Data2.fq 

Output:

	paired_MYQT${i}_Clean_Data1.fastqc.html

SourceCode:

	for i in {1..6};do
    	fastqc -o ../fastqc paired_MYQT${i}_Clean_Data1.fq 
    	fastqc -o ../fastqc paired_MYQT${i}_Clean_Data2.fq 
	done


Result:

	xxx.zip  


Note: 

	


----------------------------------------------------------------------------


Step4 Bowtie2 mapping


Input:

	GCF_000513215.1_DB11_genomic.fna
	paired_MYQT${i}_Clean_Data1.fq 
	paired_MYQT${i}_Clean_Data2.fq 
	

Output:

	MYQT${i}_Clean_Data1_fastqc.html

SourceCode (Under Checking):  

	bowtie2-build ./reference/SM_db11.fa ./reference/SM_db11.fa


	for i in {1..6};do
    	bowtie2 -p 8 -x /gsap/garage-bacterial/Users/WenChi/projects/serratiaMarcescens/data/reference/reference_RefSeq/Smarcescens -1 ./MYQT1/MYQT${i}_Clean_Data1_paried.fq.gz -2 ./MYQT1/MYQT${i}_Clean_Data2_paried.fq.gz -S ./MYQT${i}.bug.bowtie2.sam && samtools view -bS ./MYQT${i}.bug.bowtie2.sam > ./MYQT${i}.bug.bowtie2.bam  && rm ./MYQT${i}.bug.bowtie2.sam
	done


Result:

	xxx.zip


Note: 


----------

Step5 Sort bam file


Input:

	

Output:

	

SourceCode:



Result:

	


Note: 



----------

Step6 view through IGV

Input:

	GCF_000513215.1_DB11_genomic.fna
	MYQT2.bug.bowtie2.bam.onlyMapped.bam.psorted.bam
	

Output:

	

SourceCode:




Result:

	


Note: 


----------


Step7 bedtools Quantification

Input:

	GCF_000513215.1_DB11_genomic.fna
	MYQT2.bug.bowtie2.bam.onlyMapped.bam.psorted.bam
	

Output:

	

SourceCode(Under Checking):
	
	for i in {1..6};do 
		express --output-dir ./express/MYQT${i}_express ./reference/SM_db11.fa ./bam/MYQT${i}.bowtie2.bam &  
	done




Result:

	


Note: 


----------


Step8 edgeR Differential express gene 

Input:

	

Output:

	

SourceCode(Under Checking):
	
	

Result:
	


Note: 

----------


Step

Input:

	

Output:

	

SourceCode(Under Checking):
	
	

Result:
	


Note: