
########################################
######################################## for first_batch
########################################
## step1: factqc
###### qc
cd /mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure2/1.data/rawData/First_batch
nohup fastqc -o /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/2.rawfastqc -t 18 -q ./*.fastq.gz > /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/2.rawfastqc/1.fastqc.txt &
cd /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/2.rawfastqc
multiqc ./*fastqc.zip

## step2: trim R1
trimdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/3.trim
rawdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure2/1.data/rawData/First_batch
decompressdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/1.decompress
trimmomatic=/mnt/share/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimc=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R1/cutadapter
trimt=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R1/trimmomatic
trimt2=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R1/trimmomatic2
#mkdir FastQC
#mkdir FilterData
for i in WM-SGWTS-MCF7-CPT-S2_S3_L001 WM-SGWTS-MCF7-ETO-S3_S2_L001 WM-SGWTS-MCF7-IR-S3_S1_L001 WM_SGWTS_293T_S_S1_L001 WM_SGWTS_MCF7_IR_S1_L001;
do 
	#fastqc ./rawData/${i}_R1_001.fastq.gz -o ./FastQC
	gzip -d -c ${rawdir}/${i}_R1_001.fastq.gz > ${decompressdir}/${i}_R1_001.fastq
	echo "cut 5' adapter AAAATCTCTAGCA"
  	cutadapt -g ^tccggactgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagca  -e 0.2 -m 20  -o ${trimdir}/${i}_R1.filter0.1.fastq ${decompressdir}/${i}_R1_001.fastq
  	fastqc ${trimdir}/${i}_R1.filter0.3.fastq -o ${trimc}
	echo "cut 3' adapter and do quality control"
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R1.filter0.3.fastq ${trimdir}/${i}_R1.filter2.fastq ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/illumina_adapter.fa:2:30:7 CROP:150 SLIDINGWINDOW:4:20  MINLEN:20
        fastqc ${trimdir}/${i}_R1.filter2.fastq -o ${trimt}
	#rm ${trimdir}/${i}_R1.filter.fastq
	echo "cut bait sequences in 3'end"
        java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R1.filter2.fastq ${trimdir}/${i}_R1.filter2new.fastq  ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/Read1_bait_adapter.txt:2:30:7 CROP:150 MINLEN:20
         fastqc ${trimdir}/${i}_R1.filter2new.fastq -o ${trimt2}
done

## step3: trim R2
trimdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/3.trim
rawdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure2/1.data/rawData/First_batch
decompressdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/1.decompress
trimmomatic=/mnt/share/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimc=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R2/cutadapter
trimc2=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R2/cutadapter2
trimt=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R2/trimmomatic
trimt2=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/4.trimfasctqc/R2/trimmomatic2
#mkdir FastQC
#mkdir FilterData
for i in WM-SGWTS-MCF7-CPT-S2_S3_L001 WM-SGWTS-MCF7-ETO-S3_S2_L001 WM-SGWTS-MCF7-IR-S3_S1_L001;
do 
	#fastqc ./rawData/${i}_R1_001.fastq.gz -o ./FastQC
	gzip -d -c ${rawdir}/${i}_R2_001.fastq.gz > ${decompressdir}/${i}_R2_001.fastq
	cutadapt --trimmed-only -g ^CTGCGGTGATG -e 0.2 -m 20  -o ${trimdir}/${i}_R2.filter.fastq ${decompressdir}/${i}_R2_001.fastq
        fastqc ${trimdir}/${i}_R2.filter.fastq -o ${trimc}
	#cut 5' adapter CCGGTG
	cutadapt -g CCGGTG -n 6 -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter1.fastq ${trimdir}/${i}_R2.filter.fastq
		fastqc ${trimdir}/${i}_R2.filter1.fastq -o ${trimc2}
	cutadapt -g ACTCGATCTC -n 6 -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter1.5.fastq ${trimdir}/${i}_R2.filter1.fastq
	  fastqc ${trimdir}/${i}_R2.filter1.5.fastq -o ${trimc2}
 	#do quality control
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R2.filter1.5.fastq ${trimdir}/${i}_R2.filter2.fastq CROP:150 SLIDINGWINDOW:4:20 MINLEN:20
	fastqc ${trimdir}/${i}_R2.filter2.fastq -o ${trimt}
	#rm ${i}_R2.filter1.fastq
	#rm ${i}_R2.filter.fastq
	echo "cut bait sequences in 3'end"
        java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R2.filter2.fastq ${trimdir}/${i}_R2.filter2new.fastq  ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/Read2_bait_adapter.txt:2:30:7 CROP:150 MINLEN:20
   fastqc ${trimdir}/${i}_R2.filter2new.fastq -o ${trimt2}
done

## step4: find translocation and bait
trimdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/3.trim
refstargenome=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/genomeref/IndexGencodev149
starres=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/5.star

for fileName in WM-SGWTS-MCF7-CPT-S2_S3_L001 WM-SGWTS-MCF7-ETO-S3_S2_L001 WM-SGWTS-MCF7-IR-S3_S1_L001 WM_SGWTS_293T_S_S1_L001 WM_SGWTS_MCF7_IR_S1_L001;do
	echo ${fileName}
	#only keep reads in both pair
	grep '^@M' ${trimdir}/${fileName}_R1.filter2new.fastq > R1.name
	grep '^@M' ${trimdir}/${fileName}_R2.filter2new.fastq > R2.name
	cat R1.name R2.name |sed 's/ /\t/g'|cut -f 1|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}{if($1==2){print}}' >overlap.name
	wc -l overlap.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' overlap.name R1.name > R1_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' overlap.name R2.name > R2_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir}/${fileName}_R1.filter2.fastq R1_filter.name > ${trimdir}/${fileName}_R1_001.filter3.fastq
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir}/${fileName}_R2.filter2.fastq R2_filter.name > ${trimdir}/${fileName}_R2_001.filter3.fastq
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	echo "only consider the best chimeric and the ChimScoreDropMax<=20 and unique supported reads"
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}Chimeric.out.junction|cut -f 10|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g' > ${starres}/${fileName}Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}Chim_num - > ${starres}/${fileName}Chimeric.out.junction.filter
	echo "only consider read1"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	echo "only consider read2"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	echo "obtain interacted read pairs for unique mapped reads"
        samtools view -h -b -q 255 ${starres}/${fileName}_read1Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read1_unique.bam
        samtools view -h -b -q 255 ${starres}/${fileName}_read2Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read2_unique.bam
        samtools index ${starres}/${fileName}_read1_unique.bam
				samtools index ${starres}/${fileName}_read2_unique.bam     
        echo "combine reads and do filteration"
        samtools view ${starres}/${fileName}_read2_unique.bam >${starres}/${fileName}_read2_unique.sam
        samtools view ${starres}/${fileName}_read1_unique.bam >${starres}/${fileName}_read1_unique.sam
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$0}NR>FNR{if(A[$1]!=""){print $0,A[$1]}}' ${starres}/${fileName}_read2_unique.sam  ${starres}/${fileName}_read1_unique.sam |cut -f 1-6,21-25|awk 'BEGIN{FS=OFS="\t"}{if(($2!=$7) && ((($4-$9)<50 && ($4-$9)>=0) || (($4-$9)>(-50) && ($4-$9)<=0))){print}}' > ${starres}/temp
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$0}NR>FNR{if(A[$1]!=""){print $0,A[$1]}}' ${starres}/${fileName}_read2_unique.sam  ${starres}/${fileName}_read1_unique.sam |cut -f 1-6,21-25|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$1]==""){print}}' ${starres}/temp - > ${starres}/${fileName}_combine_unique.sam
        echo "obtain part1 and part3 reads"
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$2"\t"$3"\t"$4"\t"$6}NR>FNR{if(A[$10]!="" && $16-$18<=20){print $0,A[$10]}}' ${starres}/${fileName}_read2_unique.sam ${starres}/${fileName}_read1Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{dist="";if($4==$22){dist=$5-$23};if(dist!="" && ((dist>=0 && dist<=1000)||(dist<0 && dist>(-1000)))){print}}' > ${starres}/${fileName}_part1_temp
        cut -f 10 ${starres}/${fileName}_part1_temp|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1==1){A[$2]=$1}}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1_temp > ${starres}/${fileName}_part1
        echo "part2 and part3 reads"
				awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$2"\t"$3"\t"$4"\t"$6}NR>FNR{if(A[$10]!="" && $16-$18<=20){print $0,A[$10]}}' ${starres}/${fileName}_read1_unique.sam ${starres}/${fileName}_read2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{dist="";if($4==$22){dist=$5-$23};if(dist!="" && ((dist>=0 && dist<=1000)||(dist<0 && dist>(-1000)))){print}}' > ${starres}/${fileName}_part2_temp
        cut -f 10 ${starres}/${fileName}_part2_temp|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1==1){A[$2]=$1}}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part2_temp > ${starres}/${fileName}_part2
        echo "combine junction reads together"
        cat ${starres}/${fileName}_part1 ${starres}/${fileName}_part2 ${starres}/${fileName}Chimeric.out.junction.filter|cut -f 10|sort|uniq|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]="part1"}NR>FNR{if(A[$1]){print $0,A[$1]}else{print $0,"nonPart1"}}' ${starres}/${fileName}_part1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]="part2"}NR>FNR{if(A[$1]){print $0,A[$1]}else{print $0,"nonPart2"}}' ${starres}/${fileName}_part2 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($10!=""){A[$10]="PE"}}NR>FNR{if(A[$1]){print $0,A[$1]}else if($1!=""){print $0,"nonPE"}}' ${starres}/${fileName}Chimeric.out.junction.filter - > ${starres}/${fileName}_allJunctionReads
        grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part2 > part1
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1 > part2
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($2=="nonPart1" && $3=="nonPart2" && $4=="PE"){A[$1]="yes"}else if($2=="part1" && $3=="part2" && $4=="PE"){A[$1]="yes"}}NR>FNR{if(A[$10]!=""){print}}' ${starres}/${fileName}_allJunctionReads ${starres}/${fileName}Chimeric.out.junction.filter >part3
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_part1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_part2 -| awk 'BEGIN{FS=OFS="\t"}{if($8==$11 && $5==$14){dist1=$6-$15;dist2=$9-$12;if(((dist1>=0 && dist1<=150) || (dist1<=0 && dist1>=(-150)))&&((dist2>=0 && dist2<=150) || (dist2<=0 && dist2>=(-150)))){print}}}'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1 > part4
        cat part1 part2 part3 part4 |cut -f 1-14 > ${starres}/${fileName}_allJunctionReads_use
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|wc -l
        wc -l part4
        #rm part1
        #rm part2
        #rm part3
        #rm part4
				grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|grep -w 'PE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]){print}}' - ${starres}/${fileName}_allJunctionReads_use |cut -f 10,12,14|grep -v 'p'|wc -l
done

########################################
######################################## for second_batch
########################################

## step1: factqc
###### qc
cd /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.rawdata
nohup fastqc -o /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.rawfastqc -t 18 -q ./*.fastq.gz > /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.rawfastqc/1.fastqc.txt &
cd /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.rawfastqc
multiqc ./*fastqc.zip

## step2: trim R1
trimdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/3.trim
rawdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure2/1.data/rawData/Second_batch
decompressdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/1.decompress
trimmomatic=/mnt/share/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimc=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R1/cutadapter
trimt=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R1/trimmomatic
trimt2=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R1/trimmomatic2
#mkdir FastQC
#mkdir FilterData
for i in CPT21_S7_L001 ETO21_S8_L001 ETO23_S9_L001 wm-CPT-S2-R2_S3_L001 wm-ETO-S3-R2_S1_L001 wm-IR-S3-R2_S2_L001 wm-IR_S6_L001;
do 
	#fastqc ./rawData/${i}_R1_001.fastq.gz -o ./FastQC
	gzip -d -c ${rawdir}/${i}_R1_001.fastq.gz > ${decompressdir}/${i}_R1_001.fastq
	echo "cut 5' adapter AAAATCTCTAGCA"
  	cutadapt -g AAAATCTCTAGCA  -e 0.1 -m 20  -o ${trimdir}/${i}_R1.filter.fastq ${decompressdir}/${i}_R1_001.fastq
  	fastqc ${trimdir}/${i}_R1.filter.fastq -o ${trimc}
	echo "cut 3' adapter and do quality control"
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R1.filter.fastq ${trimdir}/${i}_R1.filter2.fastq ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/illumina_adapter.fa:2:30:7 CROP:150 SLIDINGWINDOW:4:20  MINLEN:20
        fastqc ${trimdir}/${i}_R1.filter2.fastq -o ${trimt}
	#rm ${trimdir}/${i}_R1.filter.fastq
	echo "cut bait sequences in 3'end"
        java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R1.filter2.fastq ${trimdir}/${i}_R1.filter2new.fastq  ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/Read1_bait_adapter.txt:2:30:7 CROP:150 MINLEN:20
         fastqc ${trimdir}/${i}_R1.filter2new.fastq -o ${trimt2}
done

## step3: trim R2
trimdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/3.trim
rawdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure2/1.data/rawData/Second_batch
decompressdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/1.decompress
trimmomatic=/mnt/share/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimc=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R2/cutadapter
trimc2=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R2/cutadapter2
trimt=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R2/trimmomatic
trimt2=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/4.trimfasctqc/R2/trimmomatic2
#mkdir FastQC
#mkdir FilterData
for i in CPT21_S7_L001 ETO21_S8_L001 ETO23_S9_L001 wm-CPT-S2-R2_S3_L001 wm-ETO-S3-R2_S1_L001 wm-IR-S3-R2_S2_L001 wm-IR_S6_L001;
do 
	#fastqc ./rawData/${i}_R1_001.fastq.gz -o ./FastQC
	gzip -d -c ${rawdir}/${i}_R2_001.fastq.gz > ${decompressdir}/${i}_R2_001.fastq
	cutadapt --trimmed-only -g ^GATAGGGATAT -e 0.1 -m 20  -o ${trimdir}/${i}_R2.filter.fastq ${decompressdir}/${i}_R2_001.fastq   ## 调整了最后的A变成T，e从0.2变成0.1
    fastqc ${trimdir}/${i}_R2.filter.fastq -o ${trimc}
	#cut 5' adapter CCGGTG CCGGTGATGCGGCACTCGATCTC
	cutadapt -g CCGGTG -n 6 -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter1.fastq ${trimdir}/${i}_R2.filter.fastq
		fastqc ${trimdir}/${i}_R2.filter1.fastq -o ${trimc2}
	cutadapt -g ACTCGATCTC -n 6 -e 0.1 -m 20 -o ${trimdir}/${i}_R2.filter1.5.fastq ${trimdir}/${i}_R2.filter1.fastq
	  fastqc ${trimdir}/${i}_R2.filter1.5.fastq -o ${trimc2}
 	#do quality control
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R2.filter1.5.fastq ${trimdir}/${i}_R2.filter2.fastq CROP:150 SLIDINGWINDOW:4:20 MINLEN:20
	fastqc ${trimdir}/${i}_R2.filter2.fastq -o ${trimt}
	#rm ${i}_R2.filter1.fastq
	#rm ${i}_R2.filter.fastq
	echo "cut bait sequences in 3'end"
        java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir}/${i}_R2.filter2.fastq ${trimdir}/${i}_R2.filter2new.fastq  ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/Read2_bait_adapter.txt:2:30:7 CROP:150 MINLEN:20
   fastqc ${trimdir}/${i}_R2.filter2new.fastq -o ${trimt2}
done

## step4: find translocation and bait
trimdir=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/3.trim
refstargenome=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/genomeref/IndexGencodev149
starres=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/5.star

for fileName in CPT21_S7_L001 ETO21_S8_L001 ETO23_S9_L001 wm-CPT-S2-R2_S3_L001 wm-ETO-S3-R2_S1_L001 wm-IR-S3-R2_S2_L001 wm-IR_S6_L001;do
	echo ${fileName}
	#only keep reads in both pair
	grep '^@M' ${trimdir}/${fileName}_R1.filter2new.fastq > R1.name
	grep '^@M' ${trimdir}/${fileName}_R2.filter2new.fastq > R2.name
	cat R1.name R2.name |sed 's/ /\t/g'|cut -f 1|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}{if($1==2){print}}' >overlap.name
	wc -l overlap.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' overlap.name R1.name > R1_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' overlap.name R2.name > R2_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir}/${fileName}_R1.fpailter2.fastq R1_filter.name > ${trimdir}/${fileName}_R1_001.filter3.fastq
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir}/${fileName}_R2.filter2.fastq R2_filter.name > ${trimdir}/${fileName}_R2_001.filter3.fastq
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	echo "only consider the best chimeric and the ChimScoreDropMax<=20 and unique supported reads"
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}Chimeric.out.junction|cut -f 10|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g' > ${starres}/${fileName}Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}Chim_num - > ${starres}/${fileName}Chimeric.out.junction.filter
	echo "only consider read1"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	echo "only consider read2"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	echo "obtain interacted read pairs for unique mapped reads"
        samtools view -h -b -q 255 ${starres}/${fileName}_read1Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read1_unique.bam
        samtools view -h -b -q 255 ${starres}/${fileName}_read2Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read2_unique.bam
        samtools index ${starres}/${fileName}_read1_unique.bam
				samtools index ${starres}/${fileName}_read2_unique.bam     
        echo "combine reads and do filteration"
        samtools view ${starres}/${fileName}_read2_unique.bam >${starres}/${fileName}_read2_unique.sam
        samtools view ${starres}/${fileName}_read1_unique.bam >${starres}/${fileName}_read1_unique.sam
        #### 输出匹配的reads
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$0}NR>FNR{if(A[$1]!=""){print $0,A[$1]}}' ${starres}/${fileName}_read2_unique.sam  ${starres}/${fileName}_read1_unique.sam |cut -f 1-6,21-25|awk 'BEGIN{FS=OFS="\t"}{if(($2!=$7) && ((($4-$9)<50 && ($4-$9)>=0) || (($4-$9)>(-50) && ($4-$9)<=0))){print}}' > ${starres}/temp
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$0}NR>FNR{if(A[$1]!=""){print $0,A[$1]}}' ${starres}/${fileName}_read2_unique.sam  ${starres}/${fileName}_read1_unique.sam |cut -f 1-6,21-25|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$1]==""){print}}' ${starres}/temp - > ${starres}/${fileName}_combine_unique.sam
        echo "obtain part1 and part3 reads"
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$2"\t"$3"\t"$4"\t"$6}NR>FNR{if(A[$10]!="" && $16-$18<=20){print $0,A[$10]}}' ${starres}/${fileName}_read2_unique.sam ${starres}/${fileName}_read1Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{dist="";if($4==$22){dist=$5-$23};if(dist!="" && ((dist>=0 && dist<=1000)||(dist<0 && dist>(-1000)))){print}}' > ${starres}/${fileName}_part1_temp
        cut -f 10 ${starres}/${fileName}_part1_temp|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1==1){A[$2]=$1}}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1_temp > ${starres}/${fileName}_part1
        echo "part2 and part3 reads"
				awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$2"\t"$3"\t"$4"\t"$6}NR>FNR{if(A[$10]!="" && $16-$18<=20){print $0,A[$10]}}' ${starres}/${fileName}_read1_unique.sam ${starres}/${fileName}_read2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{dist="";if($4==$22){dist=$5-$23};if(dist!="" && ((dist>=0 && dist<=1000)||(dist<0 && dist>(-1000)))){print}}' > ${starres}/${fileName}_part2_temp
        cut -f 10 ${starres}/${fileName}_part2_temp|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1==1){A[$2]=$1}}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part2_temp > ${starres}/${fileName}_part2
        echo "combine junction reads together"
        cat ${starres}/${fileName}_part1 ${starres}/${fileName}_part2 ${starres}/${fileName}Chimeric.out.junction.filter|cut -f 10|sort|uniq|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]="part1"}NR>FNR{if(A[$1]){print $0,A[$1]}else{print $0,"nonPart1"}}' ${starres}/${fileName}_part1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]="part2"}NR>FNR{if(A[$1]){print $0,A[$1]}else{print $0,"nonPart2"}}' ${starres}/${fileName}_part2 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($10!=""){A[$10]="PE"}}NR>FNR{if(A[$1]){print $0,A[$1]}else if($1!=""){print $0,"nonPE"}}' ${starres}/${fileName}Chimeric.out.junction.filter - > ${starres}/${fileName}_allJunctionReads
        grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part2 > part1
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1 > part2
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($2=="nonPart1" && $3=="nonPart2" && $4=="PE"){A[$1]="yes"}else if($2=="part1" && $3=="part2" && $4=="PE"){A[$1]="yes"}}NR>FNR{if(A[$10]!=""){print}}' ${starres}/${fileName}_allJunctionReads ${starres}/${fileName}Chimeric.out.junction.filter >part3
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_part1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_part2 -| awk 'BEGIN{FS=OFS="\t"}{if($8==$11 && $5==$14){dist1=$6-$15;dist2=$9-$12;if(((dist1>=0 && dist1<=150) || (dist1<=0 && dist1>=(-150)))&&((dist2>=0 && dist2<=150) || (dist2<=0 && dist2>=(-150)))){print}}}'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1 > part4
        cat part1 part2 part3 part4 |cut -f 1-14 > ${starres}/${fileName}_allJunctionReads_use
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|wc -l
        wc -l part4
        #rm part1
        #rm part2
        #rm part3
        #rm part4
				grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|grep -w 'PE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]){print}}' - ${starres}/${fileName}_allJunctionReads_use |cut -f 10,12,14|grep -v 'p'|wc -l
done

########################################
######################################## combine two batches and merge baits
########################################

R CMD BATCH ./1.combine1.r
## outputs: merged bait sites, translocation sites after removing duplications

########################################
######################################## for MCF-7 without drug treatment two replicates
########################################

## step1: factqc
cd /mnt1/2.NAS2024/share/data/WM_022024
nohup fastqc -o /mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/1.fastqc -t 18 -q ./*DOX*.fastq.gz > /mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/1.fastqc.txt &

## step2: trim R1
rawdir=/mnt1/2.NAS2024/share/data/WM_022024
decompressdir=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/decompress
trimmomatic=/mnt/share/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimdir1=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s1
trimdir2=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s2
trimdir3=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s3
trimdir4=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s4
#mkdir FastQC
#mkdir FilterData
for i in wm-DOX1_S4_L001 wm-DOX2_S5_L001;
do 
	#fastqc ./rawData/${i}_R1_001.fastq.gz -o ./FastQC
	gzip -d -c ${rawdir}/${i}_R1_001.fastq.gz > ${decompressdir}/${i}_R1_001.fastq
	echo "cut 5' adapter AAAATCTCTAGCA"
  cutadapt -g AAAATCTCTAGCA  -e 0.1 -m 20  -o ${trimdir1}/${i}_R1.filter1.fastq ${decompressdir}/${i}_R1_001.fastq   ## 几乎没有这个adapter
  #fastqc ${trimdir}/${i}_R1.filter.fastq -o ${trimc}
  echo "cut polyA"
  cutadapt -a AAAAAAAA -m 20 -o ${trimdir2}/${i}_R1.filter2.fastq ${trimdir1}/${i}_R1.filter1.fastq    ## 57%
	echo "cut 3' adapter and do quality control"
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir2}/${i}_R1.filter2.fastq ${trimdir3}/${i}_R1.filter3.fastq ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/illumina_adapter.fa:2:30:7 CROP:150 SLIDINGWINDOW:4:20  MINLEN:20
  fastqc ${trimdir3}/${i}_R1.filter3.fastq -o ${trimdir3}
	echo "cut bait sequences in 3'end"
  java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir3}/${i}_R1.filter3.fastq ${trimdir4}/${i}_R1.filter4.fastq  ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/Read1_bait_adapter.txt:2:30:7 CROP:150 MINLEN:20   ## 52%
  fastqc ${trimdir4}/${i}_R1.filter4.fastq -o ${trimdir4}
done

## step3: trim R2
rawdir=/mnt1/2.NAS2024/share/data/WM_022024
decompressdir=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/decompress
trimmomatic=/mnt/share/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
trimdir1=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s1
trimdir2=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s2
trimdir3=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s3
trimdir4=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s4
#mkdir FastQC
#mkdir FilterData
for i in wm-DOX1_S4_L001 wm-DOX2_S5_L001;
do 
	gzip -d -c ${rawdir}/${i}_R2_001.fastq.gz > ${decompressdir}/${i}_R2_001.fastq
	cutadapt --trimmed-only -g ^GATAGGGATAA -e 0.2 -m 20  -o ${trimdir1}/${i}_R2.filter.fastq ${decompressdir}/${i}_R2_001.fastq
  fastqc ${trimdir1}/${i}_R2.filter.fastq -o ${trimdir1}
	#cut 5' adapter CCGGTG
	cutadapt -g CCGGTG -n 6 -e 0.1 -m 20 -o ${trimdir1}/${i}_R2.filter1.fastq ${trimdir1}/${i}_R2.filter.fastq
	fastqc ${trimdir1}/${i}_R2.filter1.fastq -o ${trimdir1}
	cutadapt -g ACTCGATCTC -n 6 -e 0.1 -m 20 -o ${trimdir1}/${i}_R2.filter1.5.fastq ${trimdir1}/${i}_R2.filter1.fastq
	fastqc ${trimdir1}/${i}_R2.filter1.5.fastq -o ${trimdir1}
  ## cut3'
  cutadapt -a GATCGTCGGACT -m 20 -o ${trimdir2}/${i}_R2.filter2.fastq ${trimdir1}/${i}_R2.filter1.5.fastq
  cutadapt -a AAAAAAAA -m 20 -o ${trimdir2}/${i}_R2.filter2.5.fastq ${trimdir2}/${i}_R2.filter2.fastq
  fastqc ${trimdir2}/${i}_R2.filter2.5.fastq -o ${trimdir2}
 	#do quality control
	java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir2}/${i}_R2.filter2.5.fastq ${trimdir3}/${i}_R2.filter3.fastq CROP:150 SLIDINGWINDOW:4:20 MINLEN:20
	fastqc ${trimdir3}/${i}_R2.filter3.fastq -o ${trimdir3}
	#rm ${i}_R2.filter1.fastq
	#rm ${i}_R2.filter.fastq
	echo "cut bait sequences in 3'end"
  java -jar ${trimmomatic} SE -phred33 -threads 10 ${trimdir3}/${i}_R2.filter3.fastq ${trimdir4}/${i}_R2.filter4.fastq  ILLUMINACLIP:/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/Files/Read2_bait_adapter.txt:2:30:7 CROP:150 MINLEN:20
  fastqc ${trimdir4}/${i}_R2.filter4.fastq -o ${trimdir4}
done

####### step4: find pairs and merge bait within 1kb, remove duplications

trimdir=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim
refstargenome=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/genomeref/IndexGencodev149
starres=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/3.star
trimdir3=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s3
trimdir4=/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/2.trim/R1s4
cd ${trimdir}
for fileName in wm-DOX1_S4_L001 wm-DOX2_S5_L001;do
	echo ${fileName}
	#only keep reads in both pair
	grep '^@M' ${trimdir4}/${fileName}_R1.filter4.fastq > R1.name
	grep '^@M' ${trimdir4}/${fileName}_R2.filter4.fastq > R2.name
	cat R1.name R2.name |sed 's/ /\t/g'|cut -f 1|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}{if($1==2){print}}' >overlap.name
	wc -l overlap.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' overlap.name R1.name > R1_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]="yes"}NR>FNR{split($1,a," ");if(A[a[1]]!=""){print}}' overlap.name R2.name > R2_filter.name
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir3}/${fileName}_R1.filter3.fastq R1_filter.name > ${trimdir4}/${fileName}_R1_001.filter5.fastq
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if(NR%4==1){name=$1}else if(NR%4==2){A[name]=$1}else if(NR%4==3){B[name]=$1}else if(NR%4==0){C[name]=$1}}NR>FNR{if(A[$1]!=""){print $1"\n"A[$1]"\n"B[$1]"\n"C[$1]}}' ${trimdir3}/${fileName}_R2.filter3.fastq R2_filter.name > ${trimdir4}/${fileName}_R2_001.filter5.fastq
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir4}/${fileName}_R1_001.filter5.fastq ${trimdir4}/${fileName}_R2_001.filter5.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 --outTmpDir /mnt/share/FileTransfer/tmp/${fileName}_S12
	echo "only consider the best chimeric and the ChimScoreDropMax<=20 and unique supported reads"
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}Chimeric.out.junction|cut -f 10|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g' > ${starres}/${fileName}Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}Chim_num - > ${starres}/${fileName}Chimeric.out.junction.filter
	echo "only consider read1"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir4}/${fileName}_R1_001.filter5.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 --outTmpDir /mnt/share/FileTransfer/tmp/${fileName}_S1
	echo "only consider read2"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir4}/${fileName}_R2_001.filter5.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 --outTmpDir /mnt/share/FileTransfer/tmp/${fileName}_S2
	echo "obtain interacted read pairs for unique mapped reads"
        samtools view -h -b -q 255 ${starres}/${fileName}_read1Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read1_unique.bam
        samtools view -h -b -q 255 ${starres}/${fileName}_read2Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read2_unique.bam
        samtools index ${starres}/${fileName}_read1_unique.bam
				samtools index ${starres}/${fileName}_read2_unique.bam     
        echo "combine reads and do filteration"
        samtools view ${starres}/${fileName}_read2_unique.bam >${starres}/${fileName}_read2_unique.sam
        samtools view ${starres}/${fileName}_read1_unique.bam >${starres}/${fileName}_read1_unique.sam
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$0}NR>FNR{if(A[$1]!=""){print $0,A[$1]}}' ${starres}/${fileName}_read2_unique.sam  ${starres}/${fileName}_read1_unique.sam |cut -f 1-6,21-25|awk 'BEGIN{FS=OFS="\t"}{if(($2!=$7) && ((($4-$9)<50 && ($4-$9)>=0) || (($4-$9)>(-50) && ($4-$9)<=0))){print}}' > ${starres}/temp
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$0}NR>FNR{if(A[$1]!=""){print $0,A[$1]}}' ${starres}/${fileName}_read2_unique.sam  ${starres}/${fileName}_read1_unique.sam |cut -f 1-6,21-25|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$1]==""){print}}' ${starres}/temp - > ${starres}/${fileName}_combine_unique.sam
        echo "obtain part1 and part3 reads"
        awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$2"\t"$3"\t"$4"\t"$6}NR>FNR{if(A[$10]!="" && $16-$18<=20){print $0,A[$10]}}' ${starres}/${fileName}_read2_unique.sam ${starres}/${fileName}_read1Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{dist="";if($4==$22){dist=$5-$23};if(dist!="" && ((dist>=0 && dist<=1000)||(dist<0 && dist>(-1000)))){print}}' > ${starres}/${fileName}_part1_temp
        cut -f 10 ${starres}/${fileName}_part1_temp|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1==1){A[$2]=$1}}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1_temp > ${starres}/${fileName}_part1
        echo "part2 and part3 reads"
				awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=$2"\t"$3"\t"$4"\t"$6}NR>FNR{if(A[$10]!="" && $16-$18<=20){print $0,A[$10]}}' ${starres}/${fileName}_read1_unique.sam ${starres}/${fileName}_read2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}{dist="";if($4==$22){dist=$5-$23};if(dist!="" && ((dist>=0 && dist<=1000)||(dist<0 && dist>(-1000)))){print}}' > ${starres}/${fileName}_part2_temp
        cut -f 10 ${starres}/${fileName}_part2_temp|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($1==1){A[$2]=$1}}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part2_temp > ${starres}/${fileName}_part2
        echo "combine junction reads together"
        cat ${starres}/${fileName}_part1 ${starres}/${fileName}_part2 ${starres}/${fileName}Chimeric.out.junction.filter|cut -f 10|sort|uniq|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]="part1"}NR>FNR{if(A[$1]){print $0,A[$1]}else{print $0,"nonPart1"}}' ${starres}/${fileName}_part1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]="part2"}NR>FNR{if(A[$1]){print $0,A[$1]}else{print $0,"nonPart2"}}' ${starres}/${fileName}_part2 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($10!=""){A[$10]="PE"}}NR>FNR{if(A[$1]){print $0,A[$1]}else if($1!=""){print $0,"nonPE"}}' ${starres}/${fileName}Chimeric.out.junction.filter - > ${starres}/${fileName}_allJunctionReads
        #only keep reads in both pair
        grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part2 |awk 'BEGIN{FS=OFS="\t"}{print $4,$23,$6,$1,$11,$3,$10}'> Part2
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1 |awk 'BEGIN{FS=OFS="\t"}{print $1,$11,$3,$4,$23,$6,$10}'> Part1
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($2=="nonPart1" && $3=="nonPart2" && $4=="PE"){A[$1]="yes"}else if($2=="part1" && $3=="part2" && $4=="PE"){A[$1]="yes"}}NR>FNR{if(A[$10]!=""){print}}' ${starres}/${fileName}_allJunctionReads ${starres}/${fileName}Chimeric.out.junction.filter|awk 'BEGIN{FS=OFS="\t"}{print $1,$11,$3,$4,$13,$6,$10}' >PE_1
        grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_part1 -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_part2 -| awk 'BEGIN{FS=OFS="\t"}{if($8==$11 && $5==$14){dist1=$6-$15;dist2=$9-$12;if(((dist1>=0 && dist1<=150) || (dist1<=0 && dist1>=(-150)))&&((dist2>=0 && dist2<=150) || (dist2<=0 && dist2>=(-150)))){print}}}'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$10]!=""){print}}' - ${starres}/${fileName}_part1 |awk 'BEGIN{FS=OFS="\t"}{print $1,$11,$3,$4,$13,$6,$10}'> PE_2
	cat Part1 Part2 PE_1 PE_2 > ${starres}/${fileName}_baitAndTrans
	rm Part1
	rm Part2
	rm PE_1
	rm PE_2
done

########################################
######################################## combine two replicates for MCF-7 without drug treatment
########################################

R CMD BATCH ./2.combine1.r
## outputs: merged bait sites, translocation sites after removing duplications


