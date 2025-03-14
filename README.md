# High Multiplex Translocation Sequencing (HMTS)
This is the pipeline of downstream analysis for HMTS data.

## Part I Introduction
### i. Library construction
Here stands an schematic diagram illustrating HMTS.

![图片](https://github.com/user-attachments/assets/6bbd7a57-6b77-4713-919f-7d13d844017d)

There are four types of the translocation reads: \
(1) Short translocation & Long bait; \
(2) Long translocation & Short bait; \
(3) Short translocation & Short bait; \
(4) Long translocation & Long bait; \
Only the (4) translocation could be identified.

### ii. Workflow

![图片](https://github.com/user-attachments/assets/2f151fbf-8056-44f1-b1bc-d6ae0339bda8)

As illustrated in the figure, \
(i) yellow circles represent the steps where commands need to be entered; \
(ii) blue boxes represent the filtering criteria at each step. 

Here we obtained all bait and translocation pairs.

## Part II Codes for Analysis

### i. Raw Data Quality Check and Trimming

Trimming R1 data
```
rawdir=XXX         ## Enter the full path of the directory containing the raw fastq.gz data
decompressdir=XXX  ## Enter the full path of the directory to store the decompressed raw fastq.gz data
trimmomatic=XXX    ## Enter the full path of the trimmomatic.jar
trimdir=XXX        ## Enter the full path of the directory to store the trimmed fastq.gz data
Read1_bait_adapter=./Read1_bait_adapter.txt
for i in samplename1 samplename2;   ## Replace with the prefix of the raw fastq.gz data
do 
	fastqc ${rawdir}/${i}_R1_001.fastq.gz -o ${trimdir}
	gzip -d -c ${rawdir}/${i}_R1_001.fastq.gz > ${decompressdir}/${i}_R1_001.fastq
	echo "cut 5' adapter AAAATCTCTAGCA"
	cutadapt -g AAAATCTCTAGCA -e 0.1 -m 20 -o ${trimdir}/${i}_R1.filter1.fastq ${decompressdir}/${i}_R1_001.fastq
	fastqc ${trimdir}/${i}_R1.filter1.fastq -o ${trimdir}
	echo "cut 3' adapter and do quality control"
	cutadapt -g A{10} -g G{10} -a A{10} -a G{10} -a file:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/illumina_adapter.fa -j 10 --times 20 -e 0.1 -m 20 -o ${trimdir}/${i}_R1.filter2.fastq ${trimdir}/${i}_R1.filter1.fastq 
	fastqc ${trimdir}/${i}_R1.filter2.fastq -o ${trimdir}
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R1.filter2.fastq ${trimdir}/${i}_R1.filter2.5.fastq CROP:150 SLIDINGWINDOW:4:20 MINLEN:20
	fastqc ${trimdir}/${i}_R1.filter2.5.fastq -o ${trimdir}
	echo "cut bait sequences in 3'end"
	cutadapt -a file:${Read1_bait_adapter} -j 10 -e 0.1 -m 20 -O 6 -o ${trimdir}/${i}_R1.filter2new.fastq ${trimdir}/${i}_R1.filter2.5.fastq
	fastqc ${trimdir}/${i}_R1.filter2new.fastq -o ${trimdir}
done
```

Trimming R2 data
```
rawdir=XXX         ## Enter the full path of the directory containing the raw fastq.gz data
decompressdir=XXX  ## Enter the full path of the directory to store the decompressed raw fastq.gz data
trimmomatic=XXX    ## Enter the full path of the trimmomatic.jar
trimdir=XXX        ## Enter the full path of the directory to store the trimmed fastq.gz data
Read2_bait_adapter=./Read2_bait_adapter.txt
for i in samplename1 samplename2;   ## Replace with the prefix of the raw fastq.gz data
do
	fastqc ${rawdir}/${i}_R2_001.fastq.gz -o ${trimdir}
	gzip -d -c ${rawdir}/${i}_R2_001.fastq.gz > ${decompressdir}/${i}_R2_001.fastq
	echo "cut 5' adapter"
	cutadapt --trimmed-only -g ^GATAGGGATAA -e 0.2 -m 20  -o ${trimdir}/${i}_R2.filter.fastq ${decompressdir}/${i}_R2_001.fastq
	fastqc ${trimdir}/${i}_R2.filter.fastq -o ${trimdir}
	cutadapt -g CCGGTG -n 6 -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter1.fastq ${trimdir}/${i}_R2.filter.fastq
	fastqc ${trimdir}/${i}_R2.filter1.fastq -o ${trimdir}
	cutadapt -g ACTCGATCTC -n 6 -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter1.5.fastq ${trimdir}/${i}_R2.filter1.fastq
	fastqc ${trimdir}/${i}_R2.filter1.5.fastq -o ${trimdir}
	echo "cut 3' adapter and do quality control"
	cutadapt -g A{10} -g G{10} -a A{10} -a G{10} -a file:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/illumina_adapter.fa -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter2.fastq ${trimdir}/${i}_R2.filter1.5.fastq
	fastqc ${trimdir}/${i}_R2.filter2.fastq -o ${trimdir}
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R2.filter2.fastq ${trimdir}/${i}_R2.filter2.5.fastq CROP:150 SLIDINGWINDOW:4:20 MINLEN:20
	fastqc ${trimdir}/${i}_R2.filter2.5.fastq -o ${trimdir}
	echo "cut bait sequences in 3'end"
	cutadapt -a file:${Read2_bait_adapter} -j 10 -e 0.1 -m 20 -O 6 -o ${trimdir}/${i}_R2.filter2new.fastq ${trimdir}/${i}_R2.filter2.5.fastq
	fastqc ${trimdir}/${i}_R2.filter2new.fastq -o ${trimdir}
done
```

### ii. Get pairs of bait and translocation

```
ulimit -n 65535
trimdir=XXX         ## Enter the full path of the directory containing the trimmed fastq.gz data
refstargenome=XXX   ## Enter the full path of the STAR index
starres=XXX         ## Enter the full path of the directory to store the result from STAR
cd ${starres}

for fileName in samplename1 samplename2;do   ## Replace with the prefix of the trimmed fastq.gz data
	echo ${fileName}
	echo "only keep reads in both pair"
	grep '^@' ${trimdir}/${fileName}_R1.filter2new.fastq > ${fileName}_R1.name
	grep '^@' ${trimdir}/${fileName}_R2.filter2new.fastq > ${fileName}_R2.name
	cat ${fileName}_R1.name ${fileName}_R2.name |sed 's/ /\t/g'|cut -f 1|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}{if($1==2){print}}' > ${fileName}_overlap.name
	wc -l ${fileName}_overlap.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' ${fileName}_overlap.name ${fileName}_R1.name > ${fileName}_R1_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' ${fileName}_overlap.name ${fileName}_R2.name > ${fileName}_R2_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir}/${fileName}_R1.filter2.fastq ${fileName}_R1_filter.name > ${trimdir}/${fileName}_R1_001.filter3.fastq
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir}/${fileName}_R2.filter2.fastq ${fileName}_R2_filter.name > ${trimdir}/${fileName}_R2_001.filter3.fastq
	rm ${fileName}_overlap.name ${fileName}_R1_filter.name ${fileName}_R2_filter.name ${fileName}_R1.name ${fileName}_R2.name
	
	echo "mapping"
	echo "only consider read1"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read1 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --alignIntronMax 0 --alignMatesGapMax 0 --chimSegmentMin 10 --chimScoreDropMax 400 --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --chimMainSegmentMultNmax 10 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimScoreSeparation 10 --chimMultimapScoreRange 10 --chimOutJunctionFormat 1
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2_read1 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --alignIntronMax 0 --alignMatesGapMax 0 --chimSegmentMin 10 --chimScoreDropMax 400 --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --chimMainSegmentMultNmax 10 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimScoreSeparation 10 --chimMultimapScoreRange 10 --chimOutJunctionFormat 1
	echo "only consider read2"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read2 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --alignIntronMax 0 --alignMatesGapMax 0 --chimSegmentMin 10 --chimScoreDropMax 400 --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --chimMainSegmentMultNmax 10 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimScoreSeparation 10 --chimMultimapScoreRange 10 --chimOutJunctionFormat 1
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2_read2 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --alignIntronMax 0 --alignMatesGapMax 0 --chimSegmentMin 10 --chimScoreDropMax 400 --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 10 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --chimMainSegmentMultNmax 10 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimScoreSeparation 10 --chimMultimapScoreRange 10 --chimOutJunctionFormat 1
	
	echo "get chimeric correspongsing junction and bed"
	# get unique mapped reads from read1 and read2, HI tag is alignment record within a set of potential alignments for a single read
	samtools view -h -b -q 255 ${starres}/${fileName}2_read1Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read1_unique.bam
	samtools view -h -b -q 255 ${starres}/${fileName}2_read2Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read2_unique.bam
	bedtools bamtobed -i ${starres}/${fileName}_read1_unique.bam -tag HI -cigar > ${starres}/${fileName}_read1_unique.bed
	bedtools bamtobed -i ${starres}/${fileName}_read2_unique.bam -tag HI -cigar > ${starres}/${fileName}_read2_unique.bed
 	# calculate the % of perfect matched nucleotides relative to read length, cutoff: %>0.66
 	samtools view ${starres}/${fileName}_read1_unique.bam |cut -f 1,6,10,17|awk 'BEGIN{FS=OFS="\t"}{split($4,a,":");split(a[3],b,/[a-zA-Z]+/);sum=0;for(i=1;i<=length(b);i++){sum=sum+b[i]};if(sum/length($3)>-1){print $2,length($3),sum,sum/length($3)}}' | paste ${starres}/${fileName}_read1_unique.bed - |awk 'BEGIN{FS=OFS="\t"}{if($11>=0.66){print}}'>  ${starres}/${fileName}_read1_unique.percent.bed
 	samtools view ${starres}/${fileName}_read2_unique.bam |cut -f 1,6,10,17|awk 'BEGIN{FS=OFS="\t"}{split($4,a,":");split(a[3],b,/[a-zA-Z]+/);sum=0;for(i=1;i<=length(b);i++){sum=sum+b[i]};if(sum/length($3)>-1){print $2,length($3),sum,sum/length($3)}}' | paste ${starres}/${fileName}_read2_unique.bed - |awk 'BEGIN{FS=OFS="\t"}{if($11>=0.66){print}}'>  ${starres}/${fileName}_read2_unique.percent.bed
 	# get all chimeric reads' bed
	samtools view -h ${starres}/${fileName}_read1Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}_read1.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}_read1.chimOut.bam -tag HI -cigar > ${starres}/${fileName}_read1.chimOut.bed
	samtools view -h ${starres}/${fileName}_read2Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}_read2.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}_read2.chimOut.bam -tag HI -cigar > ${starres}/${fileName}_read2.chimOut.bed
 	# 1) if chimeric score ($18) larger than non-chimeric score ($17), then we chose chimeric 2) chimeric mapping Bps($18) should larger then 0.66 percent of totoal reads length 
 	grep '^chr' ${starres}/${fileName}2_read1Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{if(NR >1 && $17<=$18 && $18/$16>=0.66){print}}' > ${starres}/${fileName}2_read1Chimeric.out.junction.filter1
 	grep '^chr' ${starres}/${fileName}2_read2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{if(NR >1 && $17<=$18 && $18/$16>=0.66){print}}' > ${starres}/${fileName}2_read2Chimeric.out.junction.filter1

	echo "dividing reads into different condition"
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(A[$4]==""){A[$4]=$0}else{A[$4]=A[$4]"&&&"$0}}NR>FNR{if(A[$10]!=""){split(A[$10],Loci,"&&&");for(i=1;i<=length(Loci);i++){print $0,Loci[i]}}}' ${starres}/${fileName}_read2_unique.percent.bed ${starres}/${fileName}2_read1Chimeric.out.junction.filter1 > ${starres}/${fileName}_R1J_R2U
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(A[$10]==""){A[$10]=$0}else{A[$10]=A[$10]"&&&"$0}}NR>FNR{if(A[$10]!=""){split(A[$10],Loci,"&&&");for(i=1;i<=length(Loci);i++){print $0,Loci[i]}}}' ${starres}/${fileName}2_read2Chimeric.out.junction.filter1 ${starres}/${fileName}2_read1Chimeric.out.junction.filter1 > ${starres}/${fileName}_R1J_R2J
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(A[$4]==""){A[$4]=$0}else{A[$4]=A[$4]"&&&"$0}}NR>FNR{if(A[$10]!=""){split(A[$10],Loci,"&&&");for(i=1;i<=length(Loci);i++){print $0,Loci[i]}}}' ${starres}/${fileName}_read1_unique.percent.bed ${starres}/${fileName}2_read2Chimeric.out.junction.filter1 > ${starres}/${fileName}_R2J_R1U
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(A[$4]==""){A[$4]=$0}else{A[$4]=A[$4]"&&&"$0}}NR>FNR{if(A[$4]!=""){split(A[$4],Loci,"&&&");for(i=1;i<=length(Loci);i++){print $0,Loci[i]}}}' ${starres}/${fileName}_read2_unique.percent.bed ${starres}/${fileName}_read1_unique.percent.bed > ${starres}/${fileName}_R1U_R2U
	
 	echo "deal with R1J and R2U"
	#first combine the bed together for the same map loci
	awk 'BEGIN{FS=OFS="\t"}{if(A[$4"\t"$5]!=""){print $0,A[$4"\t"$5]}else{A[$4"\t"$5]=$0}}' ${starres}/${fileName}_read1.chimOut.bed > ${starres}/${fileName}_read1.chimOut_pair.bed
	awk 'BEGIN{FS=OFS="\t"}{if(A[$4"\t"$5]!=""){print $0,A[$4"\t"$5]}else{A[$4"\t"$5]=$0}}' ${starres}/${fileName}_read2.chimOut.bed > ${starres}/${fileName}_read2.chimOut_pair.bed
	#annotated the map loci bed to the junction region
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{str1=$1"\t"$2"\t"$3"\t"$6"\t"$8"\t"$9"\t"$10"\t"$13;str2=$8"\t"$9"\t"$10"\t"$13"\t"$1"\t"$2"\t"$3"\t"$6;if($6=="+" && $13=="+"){A[$1$3+1$4$8$9]=str1;A[$8$10+1$4$1$2]=str2}else if($6=="+" && $13=="-"){A[$1$3+1$4$8$10+1]=str1;A[$8$9$4$1$2]=str2}else if($6=="-" && $13=="+"){A[$1$2$4$8$9]=str1;A[$8$10+1$4$1$3+1]=str2}else if($6=="-" && $13=="-"){A[$1$2$4$8$10+1]=str1;A[$8$9$4$1$3+1]=str2}}NR>FNR{if(A[$1$2$10$4$5]!=""){print $0,A[$1$2$10$4$5]}}' ${starres}/${fileName}_read1.chimOut_pair.bed ${starres}/${fileName}_R1J_R2U > ${starres}/${fileName}_R1J_R2U_addPair
	#filter according to R2
	awk 'BEGIN{FS=OFS="\t"}{if($21==$36 && $39!=$26){if($39=="+"){if(($38>=$23 && $38-$23<=1000) || ($38<=$23 && $23-$38<=1000)){print}}else{if(($37>=$22 && $37-$22<=1000) || ($37<=$22 && $22-$37<=1000)){print}}}}' ${starres}/${fileName}_R1J_R2U_addPair > ${starres}/${fileName}_R1J_R2U_addPair_filter
	#we should keep best R1J for mutiple mapping,and if have mutiple best chimeric,we should discard this read.
	awk 'BEGIN{FS=OFS="\t"}{if(A[$10]!=""){A[$10]=$18;B[$10]=$0}else{if($18>A[$10]){A[$10]=$18;B[$10]=$0}}}END{for(i in B){print B[i]}}' ${starres}/${fileName}_R1J_R2U_addPair_filter|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$18}NR>FNR{if($18==A[$10]){print}}' - ${starres}/${fileName}_R1J_R2U_addPair_filter|awk 'BEGIN{FS=OFS="\t"}{if(A[$10]==""){A[$10]=1;B[$10]=$0}else{A[$10]=A[$10]+1}}END{for(i in A){if(A[i]==1){print B[i]}}}' > ${starres}/${fileName}_R1J_R2U_addPair_filter1

	echo "deal with R1J and R2J"
	#Only keep R1J and R2J sharing the same junction sites(+/-10bp)
	awk 'BEGIN{FS=OFS="\t"}{if($1==$24 && $3!=$26 && $4==$21 && $6!=$23 && ((($25-$2)>=0 && ($25-$2)<=10) || (($2-$25)>=0 && ($2-$25)<=10)) && ((($22-$5)>=0 && ($22-$5)<=10) || (($5-$22)>=0 && ($5-$22)<=10))){print}}' ${starres}/${fileName}_R1J_R2J > ${starres}/${fileName}_R1J_R2J_filter
	#obtain the first loci of R1 and R2.
	#the file format: R1J,R2J,R1_part1,R1_part2,R2_part1,R2_part2
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{str1=$1"\t"$2"\t"$3"\t"$6"\t"$8"\t"$9"\t"$10"\t"$13;str2=$8"\t"$9"\t"$10"\t"$13"\t"$1"\t"$2"\t"$3"\t"$6;if($6=="+" && $13=="+"){A[$1$3+1$4$8$9]=str1;A[$8$10+1$4$1$2]=str2}else if($6=="+" && $13=="-"){A[$1$3+1$4$8$10+1]=str1;A[$8$9$4$1$2]=str2}else if($6=="-" && $13=="+"){A[$1$2$4$8$9]=str1;A[$8$10+1$4$1$3+1]=str2}else if($6=="-" && $13=="-"){A[$1$2$4$8$10+1]=str1;A[$8$9$4$1$3+1]=str2}}NR>FNR{if(A[$1$2$10$4$5]!=""){print $0,A[$1$2$10$4$5]}}' ${starres}/${fileName}_read1.chimOut_pair.bed ${starres}/${fileName}_R1J_R2J_filter | awk 'BEGIN{FS=OFS="\t"}NR==FNR{str1=$1"\t"$2"\t"$3"\t"$6"\t"$8"\t"$9"\t"$10"\t"$13;str2=$8"\t"$9"\t"$10"\t"$13"\t"$1"\t"$2"\t"$3"\t"$6;if($6=="+" && $13=="+"){A[$1$3+1$4$8$9]=str1;A[$8$10+1$4$1$2]=str2}else if($6=="+" && $13=="-"){A[$1$3+1$4$8$10+1]=str1;A[$8$9$4$1$2]=str2}else if($6=="-" && $13=="+"){A[$1$2$4$8$9]=str1;A[$8$10+1$4$1$3+1]=str2}else if($6=="-" && $13=="-"){A[$1$2$4$8$10+1]=str1;A[$8$9$4$1$3+1]=str2}}NR>FNR{if(A[$21$22$30$24$25]!=""){print $0,A[$21$22$30$24$25]}}' ${starres}/${fileName}_read2.chimOut_pair.bed - > ${starres}/${fileName}_R1J_R2J_filter_addPair
	####obtain each read the best chimeric pair score
	awk 'BEGIN{FS=OFS="\t"}{print ($18/$16)+($38/$36),$0}' ${starres}/${fileName}_R1J_R2J_filter_addPair >temp.txt
	awk 'BEGIN{FS=OFS="\t"}{if(A[$11]!=""){A[$11]=$1;B[$11]=$0}else{if($1>A[$11]){A[$11]=$1;B[$11]=$0}}}END{for(i in B){print B[i]}}' temp.txt |awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$11]=$1}NR>FNR{if(A[$11]==$1){print}}' - temp.txt|awk 'BEGIN{FS=OFS="\t"}{if(A[$11]==""){A[$11]=1;B[$11]=$0}else{A[$11]=A[$11]+1}}END{for(i in A){if(A[i]==1){print B[i]}}}' |cut -f 2- > ${starres}/${fileName}_R1J_R2J_addPair_filter1
	rm temp.txt

	echo "deal with R1U and R2J"
	#annotated the map loci bed to the junction region
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{str1=$1"\t"$2"\t"$3"\t"$6"\t"$8"\t"$9"\t"$10"\t"$13;str2=$8"\t"$9"\t"$10"\t"$13"\t"$1"\t"$2"\t"$3"\t"$6;if($6=="+" && $13=="+"){A[$1$3+1$4$8$9]=str1;A[$8$10+1$4$1$2]=str2}else if($6=="+" && $13=="-"){A[$1$3+1$4$8$10+1]=str1;A[$8$9$4$1$2]=str2}else if($6=="-" && $13=="+"){A[$1$2$4$8$9]=str1;A[$8$10+1$4$1$3+1]=str2}else if($6=="-" && $13=="-"){A[$1$2$4$8$10+1]=str1;A[$8$9$4$1$3+1]=str2}}NR>FNR{if(A[$1$2$10$4$5]!=""){print $0,A[$1$2$10$4$5]}}' ${starres}/${fileName}_read2.chimOut_pair.bed ${starres}/${fileName}_R2J_R1U > ${starres}/${fileName}_R2J_R1U_addPair
	#filter according to R1
	awk 'BEGIN{FS=OFS="\t"}{if($21==$36 && $39!=$26){if($39=="+"){if(($38>=$23 && $38-$23<=1000) || ($38<=$23 && $23-$38<=1000)){print}}else{if(($37>=$22 && $37-$22<=1000) || ($37<=$22 && $22-$37<=1000)){print}}}}' ${starres}/${fileName}_R2J_R1U_addPair > ${starres}/${fileName}_R2J_R1U_addPair_filter
	#we should keep best R2J for mutiple mapping,and if have mutiple best chimeric,we should discard this read.
	awk 'BEGIN{FS=OFS="\t"}{if(A[$10]!=""){A[$10]=$18;B[$10]=$0}else{if($18>A[$10]){A[$10]=$18;B[$10]=$0}}}END{for(i in B){print B[i]}}' ${starres}/${fileName}_R2J_R1U_addPair_filter|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$18}NR>FNR{if($18==A[$10]){print}}' - ${starres}/${fileName}_R2J_R1U_addPair_filter|awk 'BEGIN{FS=OFS="\t"}{if(A[$10]==""){A[$10]=1;B[$10]=$0}else{A[$10]=A[$10]+1}}END{for(i in A){if(A[i]==1){print B[i]}}}' > ${starres}/${fileName}_R2J_R1U_addPair_filter1

	echo "deal with R1U and R2U"
	#raise ratio to 0.75
	awk 'BEGIN{FS=OFS="\t"} {if ($11 >= 0.75 && $22 >= 0.75) {print $0}}' ${starres}/${fileName}_R1U_R2U > ${starres}/${fileName}_R1U_R2U_filter1
	
	echo "combine all reads from four types"
	R1J_R2J=${starres}/${fileName}_R1J_R2J_addPair_filter1
	R1J_R2U=${starres}/${fileName}_R1J_R2U_addPair_filter1
	R2J_R1U=${starres}/${fileName}_R2J_R1U_addPair_filter1
	R1U_R2U=${starres}/${fileName}_R1U_R2U_filter1
	# extract ReadName and type
	awk '{print $10 "\t1\t0\t0\t0"}' $R1J_R2J > ${fileName}_temp_R1J_R2J.txt
	awk '{print $10 "\t0\t1\t0\t0"}' $R1J_R2U > ${fileName}_temp_R1J_R2U.txt
	awk '{print $10 "\t0\t0\t1\t0"}' $R2J_R1U > ${fileName}_temp_R2J_R1U.txt
	awk '{print $4 "\t0\t0\t0\t1"}' $R1U_R2U > ${fileName}_temp_R1U_R2U.txt
	cat ${fileName}_temp_R1J_R2J.txt ${fileName}_temp_R1J_R2U.txt ${fileName}_temp_R2J_R1U.txt ${fileName}_temp_R1U_R2U.txt > ${fileName}_combined.txt
	# assign to read2type dataframe
	awk 'BEGIN {FS=OFS="\t"; print "ReadName", "R1J_R2J", "R1J_R2U", "R2J_R1U", "R1U_R2U";} {readname = $1; r1j_r2j[readname] += $2; r1j_r2u[readname] += $3; r2j_r1u[readname] += $4; r1u_r2u[readname] += $5;} END {for (readname in r1j_r2j) {print readname, r1j_r2j[readname], r1j_r2u[readname], r2j_r1u[readname], r1u_r2u[readname];}}' ${fileName}_combined.txt > ${starres}/${fileName}_Readtypes.tsv
	rm ${fileName}_temp_R1J_R2J.txt ${fileName}_temp_R1J_R2U.txt ${fileName}_temp_R2J_R1U.txt ${fileName}_temp_R1U_R2U.txt ${fileName}_combined.txt
	
	echo "deal with R1J_R2U and R2J_R1U"
	#deal with reads do not belong to R1J_R2J, but belong to both R1J_R2U and R2J_R1U, by comparing the mapped reads of chimeric + unique bps from two types.
	awk 'BEGIN{FS=OFS="\t"}{if($2==0 && $3==1 && $4==1){print}}' ${starres}/${fileName}_Readtypes.tsv|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$18+$30}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_R1J_R2U_addPair_filter1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$18+$30}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_R2J_R1U_addPair_filter1 -| awk 'BEGIN{FS=OFS="\t"}{if($6>$5){$4=0}else if($6<$5){$3=0};print}' | cut -f 1-5 > ${fileName}_R1J_R2U_R2J_R1U.txt
	awk 'BEGIN{FS=OFS="\t"} !($2==0 && $3==1 && $4==1){print}' ${starres}/${fileName}_Readtypes.tsv | cat - ${fileName}_R1J_R2U_R2J_R1U.txt > ${starres}/${fileName}_Readtypes_revised.tsv
	rm ${fileName}_R1J_R2U_R2J_R1U.txt
	
	echo "get bait and translocation"
	read2type=${starres}/${fileName}_Readtypes_revised.tsv
	# extract from R1J_R2J directly, get 5' sites from read
	awk 'BEGIN{FS=OFS="\t"} {if ($44 == "+" && $52 == "+") {print $41, $42, $42+1, $44, $49, $50, $50+1, $52, $10} else if ($44 == "+" && $52 == "-") {print $41, $42, $42+1, $44, $49, $51-1, $51, $52, $10} else if ($44 == "-" && $52 == "-") {print $41, $43-1, $43, $44, $49, $51-1, $51, $52, $10} else if ($44 == "-" && $52 == "+") {print $41, $43-1, $43, $44, $49, $50, $50+1, $52, $10}}' $R1J_R2J > ${starres}/${fileName}_R1J_R2J.tsv
	# extract R1J_R2U reads, get 5' sites from read
	awk 'BEGIN{FS=OFS="\t"} $2 == 0 && $3 == 1 {print $1}' $read2type | awk 'BEGIN{FS=OFS="\t"} NR==FNR{A[$1]="yes"} NR>FNR {if (A[$10] == "yes") {if ($35 == "+" && $26 == "+") {print $32, $33, $33+1, $35, $21, $22, $22+1, $26, $10} else if ($35 == "+" && $26 == "-") {print $32, $33, $33+1, $35, $21, $23-1, $23, $26, $10} else if ($35 == "-" && $26 == "-") {print $32, $34-1, $34, $35, $21, $23-1, $23, $26, $10} else if ($35 == "-" && $26 == "+") {print $32, $34-1, $34, $35, $21, $22, $22+1, $26, $10}}}' - $R1J_R2U > ${starres}/${fileName}_R1J_R2U.tsv
	# extract R1U_R2J reads, get 5' sites from read
	awk 'BEGIN{FS=OFS="\t"} $2 == 0 && $4 == 1 {print $1}' $read2type |	awk 'BEGIN{FS=OFS="\t"} NR==FNR{A[$1]="yes"} NR>FNR {if (A[$10] == "yes") {if ($35 == "+" && $26 == "+") {print $21, $22, $22+1, $26, $32, $33, $33+1, $35, $10} else if ($35 == "+" && $26 == "-") {print $21, $23-1, $23, $26, $32, $33, $33+1, $35, $10} else if ($35 == "-" && $26 == "-") {print $21, $23-1, $23, $26, $32, $34-1, $34, $35, $10} else if ($35 == "-" && $26 == "+") {print $21, $22, $22+1, $26, $32, $34-1, $34, $35, $10}}}' - $R2J_R1U > ${starres}/${fileName}_R2J_R1U.tsv
	# extract R1U_R2U reads, get 5' sites from read
	awk 'BEGIN{FS=OFS="\t"} $2 == 0 && $3 == 0 && $4 == 0 && $5 == 1 {print $1}' $read2type | awk 'BEGIN{FS=OFS="\t"} NR==FNR{A[$1]="yes"} NR>FNR {if (A[$4] == "yes") {if ($6 == "+" && $17 == "+") {print $1, $2, $2+1, $6, $12, $13, $13+1, $17, $4} else if ($6 == "+" && $17 == "-") {print $1, $2, $2+1, $6, $12, $14-1, $14, $17, $4} else if ($6 == "-" && $17 == "-") {print $1, $3-1, $3, $6, $12, $14-1, $14, $17, $4} else if ($6 == "-" && $17 == "+") {print $1, $3-1, $3, $6, $12, $13, $13+1, $17, $4}}}' - $R1U_R2U > ${starres}/${fileName}_R1U_R2U.tsv
	# combine
	cat ${starres}/${fileName}_R1J_R2J.tsv ${starres}/${fileName}_R1J_R2U.tsv ${starres}/${fileName}_R2J_R1U.tsv ${starres}/${fileName}_R1U_R2U.tsv > ${starres}/${fileName}_allBaitTransReads_1bp
	rm ${starres}/${fileName}_R1J_R2J.tsv ${starres}/${fileName}_R1J_R2U.tsv ${starres}/${fileName}_R2J_R1U.tsv ${starres}/${fileName}_R1U_R2U.tsv
done
```



