
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
refstargenome=/mnt/share/Index/STAR_index/hg38_chr_149
starres=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/7.star

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
	echo "paried-end mapping"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	rm ${starres}/${fileName}2Aligned.sortedByCoord.out.bam
	#
	echo "only consider the best chimeric and the ChimScoreDropMax<=20 and unique supported reads"
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}2Chimeric.out.junction|cut -f 10|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g' > ${starres}/${fileName}Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}Chim_num - > ${starres}/${fileName}Chimeric.out.junction.filter
	samtools view -h ${starres}/${fileName}Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}.chimOut.bam -tag HI -cigar > ${starres}/${fileName}.chimOut.bed
	cat ${starres}/${fileName}Chimeric.out.junction.filter | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($12 ~ /[0-9]+M/) {
        match($12, /[0-9]+M/, m);
        $21 = m[0];
    } else {
        $21 = ".";
    }
    if ($14 ~ /[0-9]+M/) {
        match($14, /[0-9]+M/, m);
        $22 = m[0];
    } else {
        $22 = ".";
    }
    print $0;
}' - > ${starres}/${fileName}Chimeric.out.junction.filter2
	cat ${starres}/${fileName}.chimOut.bed | awk '
BEGIN { OFS="\t"; }
{
    # 分割第四列，存放到第八列和第九列
    split($4, arr, "/");
    $8 = arr[1];
    $9 = arr[2];
    
    # 提取第七列中的 M 和数字
    if ($7 ~ /[0-9]+M/) {
        match($7, /[0-9]+M/, m);
        $10 = m[0];
    } else {
        $10 = ".";
    }
    
    print $0;
}' - | sort -k8,8 -k5,5n -k9,9 | awk '
BEGIN {
    OFS="\t";
}
{
    key = $8 ":" $5;
    if (key != prev_key) {
        # 如果新分组且前面有输出，则输出上一组
        if (NR > 1) {
            if (count == 2) {
                output_line = output_line "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN";
            }
            print output_line;
        }
        # 初始化新组
        output_line = $0;
        count = 1;
    } else {
        # 如果key相同，将新行串联到当前组
        count++;
        output_line = output_line "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10;
    }
    prev_key = key;
}
END {
    # 输出最后一组的数据
    if (count == 2) {
        output_line = output_line "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN";
    }
    print output_line;
}' - | awk 'BEGIN { OFS="\t"; } 
NR==FNR { aa_map[$10] = $0; next; } 
{
    if ($8 in aa_map) {
        print $0, aa_map[$8];
    } else {
        print $0, ".";
    }
}' ${starres}/${fileName}Chimeric.out.junction.filter2 - | awk '
BEGIN {
    OFS="\t";
}
function abs(x) {
    return (x < 0) ? -x : x
}
{

    # 如果第21列是 "NN"
    if ($21 == "NN") {
        if ($1 == $31 && $11 == $34 && (index($42, $10) > 0 || index($44, $10) > 0) && (index($42, $20) > 0 || index($44, $20) > 0) && (index($7, $51) > 0 || index($17, $51) > 0) && (index($7, $52) > 0 || index($17, $52) > 0) && ($32 == $2 || $32 == $12 || $32 - 1 == $3 || $32 - 1 == $13 || $32 + 1 == $3 || $32 + 1 == $13) && ($35 == $2|| $35 == $12 || $35 - 1 == $3 || $35 - 1 == $13 || $35 + 1 == $3 || $35 + 1 == $13)) {
            # 符合条件，检查第6列的符号并进行相应操作
            if ($6 == "+") {
                $3 = $2 + 1
            } else if ($6 == "-") {
                $2 = $3 - 1
            }
            if ($16 == "+") {
                $13 = $12 + 1
            } else if ($16 == "-") {
                $12 = $13 - 1
            }
            # 输出 1-7 和 11-17 列
            print $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14, $15, $16, $17, $18
        }
    }
    # 如果第21列不是 "NN"
    else {
        if ((($1 == $31 || $1 == $34) && 
             ($11 == $31 || $11 == $34) && 
             ($21 == $31 || $21 == $34) && 
             (index($42, $10) > 0 || index($44, $10) > 0) && 
             (index($42, $20) > 0 || index($44, $20) > 0) && 
             (index($42, $30) > 0 || index($44, $30) > 0) && 
             (index($7, $51) > 0 || index($17, $51) > 0 || index($27, $51) > 0) && 
             (index($7, $52) > 0 || index($17, $52) > 0 || index($27, $52) > 0) && 
             ($32 == $2 || $32 == $12 || $32 == $22 || $32 - 1 == $3 || $32 - 1 == $13 || $32 - 1 == $23 || $32 + 1 == $3 || $32 + 1 == $13 || $32 + 1 == $23) && 
             ($35 == $2 || $35 == $12 || $35 == $22 || $35 - 1 == $3 || $35 - 1 == $13 || $35 - 1 == $23 || $35 + 1 == $3 || $35 + 1 == $13 || $35 + 1 == $23))) {
            # 计算selpos
            if ($19 == 1) {
                if ($1 == $21) {
                    if ($11 == $21) {
                        if (abs($2 - $32) > abs($12 - $32)) {
                            selpos = 11
                        } else {
                            selpos = 12
                        }
                    } else {
                        selpos = 12
                    }
                } else {
                    selpos = 11
                }
            } else if ($19 == 2) {
                if ($11 == $1) {
                    if ($21 == $1) {
                        if (abs($12 - $35) > abs($22 - $35)) {
                            selpos = 22
                        } else {
                            selpos = 23
                        }
                    } else {
                        selpos = 23
                    }
                } else {
                    selpos = 22
                }
            }
            
            # 根据selpos更新列并输出不同的列
            if (selpos == 11 || selpos == 23) {
                if ($6 == "+") {
                    $3 = $2 + 1
                } else if ($6 == "-") {
                    $2 = $3 - 1
                }
                if ($26 == "+") {
                    $23 = $22 + 1
                } else if ($26 == "-") {
                    $22 = $23 - 1
                }
                # 输出1-7 和 21-27列
                print $1, $2, $3, $4, $5, $6, $7, $21, $22, $23, $24, $25, $26, $27, $28
            } else if (selpos == 12) {
                if ($16 == "+") {
                    $13 = $12 + 1
                } else if ($16 == "-") {
                    $12 = $13 - 1
                }
                if ($26 == "+") {
                    $23 = $22 + 1
                } else if ($26 == "-") {
                    $22 = $23 - 1
                }
                # 输出11-17 和 21-27列
                print $11, $12, $13, $14, $15, $16, $17, $21, $22, $23, $24, $25, $26, $27, $28
            } else if (selpos == 22) {
                if ($6 == "+") {
                    $3 = $2 + 1
                } else if ($6 == "-") {
                    $2 = $3 - 1
                }
                if ($16 == "+") {
                    $13 = $12 + 1
                } else if ($16 == "-") {
                    $12 = $13 - 1
                }
                # 输出1-7 和 11-17列
                print $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14, $15, $16, $17, $18
            }
        }
    }
}' - > ${starres}/${fileName}Chimeric.junction.5pend.txt
	echo "only consider read1"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	echo "only consider read2"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	echo "obtain interacted read pairs for unique mapped reads"
	samtools view -h -b -q 255 ${starres}/${fileName}2_read1Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read1_unique.bam
 	samtools view -h -b -q 255 ${starres}/${fileName}2_read2Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read2_unique.bam
 	bedtools bamtobed -i ${starres}/${fileName}_read1_unique.bam -tag HI -cigar > ${starres}/${fileName}_read1_unique.bed
 	bedtools bamtobed -i ${starres}/${fileName}_read2_unique.bam -tag HI -cigar > ${starres}/${fileName}_read2_unique.bed
 	echo "combine reads and do filteration for read1"
	awk 'BEGIN{FS=OFS="\t"} !/^#/ && $18==$19 && $16-$18<=20 {print}' ${starres}/${fileName}2_read1Chimeric.out.junction | cut -f 10 | sort | uniq -c | sed 's/^[t]*//g' | sed 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"} {n=0; for (i=1; i<=NF; i++) if ($i != "") {if (n > 0) printf "\t"; printf "%s", $i; n++} print ""}' > ${starres}/${fileName}_read1Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}' ${starres}/${fileName}2_read1Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}_read1Chim_num - > ${starres}/${fileName}_read1Chimeric.out.junction.filter
	samtools view -h ${starres}/${fileName}_read1Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}_read1.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}_read1.chimOut.bam -tag HI -cigar > ${starres}/${fileName}_read1.chimOut.bed      
	cat ${starres}/${fileName}_read1Chimeric.out.junction.filter | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($12 ~ /[0-9]+M/) {
        match($12, /[0-9]+M/, m);
        $21 = m[0];
    } else {
        $21 = ".";
    }
    if ($14 ~ /[0-9]+M/) {
        match($14, /[0-9]+M/, m);
        $22 = m[0];
    } else {
        $22 = ".";
    }
    print $0;
}' - > ${starres}/${fileName}_read1Chimeric.out.junction.filter2
	cat ${starres}/${fileName}_read1.chimOut.bed | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($7 ~ /[0-9]+M/) {
        match($7, /[0-9]+M/, m);
        $8 = m[0];
    } else {
        $8 = ".";
    }
    print $0;
}' - | sort -k4,4 -k5,5n | awk 'BEGIN {FS=OFS="\t"}
    {
        key = $4 ":" $5;
        if (key in lines) {
            lines[key] = lines[key] OFS $0;
        } else {
            lines[key] = $0;
        }
    }
    END {
        for (key in lines) {
            print lines[key];
        }
    }' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第10列是否匹配第一个文件的第四列
        if ($10 in file1) {
            for (i = 1; i <= count[$10]; i++) {
                print file1[$10][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read1Chimeric.out.junction.filter2 | awk '
BEGIN {
    OFS="\t";
}
{
	if ((($17 == $1 || $17 == $9)) && (($20 == $1 || $20 == $9)) && ((index($28, $8) > 0 || index($30, $8) > 0)) && ((index($28, $16) > 0 || index($30, $16) > 0)) && ($18 == $2 || $18 == $3 || $18 == $10 || $18 == $11 || $18 + 1 == $2 || $18 + 1 == $3 || $18 + 1 == $10 || $18 + 1 == $11 || $18 - 1 == $2 || $18 - 1 == $3 || $18 - 1 == $10 || $18 - 1 == $11) && ($21 == $2 || $21 == $3 || $21 == $10 || $21 == $11 || $21 + 1 == $2 || $21 + 1 == $3 || $21 + 1 == $10 || $21 + 1 == $11 || $21 - 1 == $2 || $21 - 1 == $3 || $21 - 1 == $10 || $21 - 1 == $11)) {
		print $0
	}
}' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第四列是否匹配第一个文件的第四列
        if ($4 in file1) {
            # 如果匹配，则为第二个文件的每一行打印所有第一个文件中对应第四列的行
            for (i = 1; i <= count[$4]; i++) {
                print file1[$4][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read2_unique.bed | awk '
BEGIN {
    OFS="\t";
}
function abs(x) {
    return (x < 0) ? -x : x
}
{
		if((($39 == $1 && $39 != $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) < abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $3 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $2 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $3 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $2 + 2 && abs($40 - $21) < 1000))))) {
			if ($14 == "+") {
				$11 = $10 + 1
			} else if ($14 == "-") {
				$10 = $11 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $9, $10, $11, $12, $13, $14, $15, $39, $40, $41, $42, $43, $44, $45, $42
		} else if((($39 != $1 && $39 == $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) > abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $11 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $10 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $11 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $10 + 2 && abs($40 - $21) < 1000)))))
			{
			if ($6 == "+") {
				$3 = $2 + 1
			} else if ($6 == "-") {
				$2 = $3 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $1, $2, $3, $4, $5, $6, $7, $39, $40, $41, $42, $43, $44, $45, $42
		}
}' - > ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt
	echo "combine reads and do filteration for read2"
	awk 'BEGIN{FS=OFS="\t"} !/^#/ && $18==$19 && $16-$18<=20 {print}' ${starres}/${fileName}2_read2Chimeric.out.junction | cut -f 10 | sort | uniq -c | sed 's/^[t]*//g' | sed 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"} {n=0; for (i=1; i<=NF; i++) if ($i != "") {if (n > 0) printf "\t"; printf "%s", $i; n++} print ""}' > ${starres}/${fileName}_read2Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}' ${starres}/${fileName}2_read2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}_read2Chim_num - > ${starres}/${fileName}_read2Chimeric.out.junction.filter
	samtools view -h ${starres}/${fileName}_read2Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}_read2.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}_read2.chimOut.bam -tag HI -cigar > ${starres}/${fileName}_read2.chimOut.bed
	cat ${starres}/${fileName}_read2Chimeric.out.junction.filter | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($12 ~ /[0-9]+M/) {
        match($12, /[0-9]+M/, m);
        $21 = m[0];
    } else {
        $21 = ".";
    }
    if ($14 ~ /[0-9]+M/) {
        match($14, /[0-9]+M/, m);
        $22 = m[0];
    } else {
        $22 = ".";
    }
    print $0;
}' - > ${starres}/${fileName}_read2Chimeric.out.junction.filter2
	cat ${starres}/${fileName}_read2.chimOut.bed | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($7 ~ /[0-9]+M/) {
        match($7, /[0-9]+M/, m);
        $8 = m[0];
    } else {
        $8 = ".";
    }
    print $0;
}' - | sort -k4,4 -k5,5n | awk 'BEGIN {FS=OFS="\t"}
    {
        key = $4 ":" $5;
        if (key in lines) {
            lines[key] = lines[key] OFS $0;
        } else {
            lines[key] = $0;
        }
    }
    END {
        for (key in lines) {
            print lines[key];
        }
    }' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第10列是否匹配第一个文件的第四列
        if ($10 in file1) {
            for (i = 1; i <= count[$10]; i++) {
                print file1[$10][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read2Chimeric.out.junction.filter2 | awk '
BEGIN {
    OFS="\t";
}
{
	if ((($17 == $1 || $17 == $9)) && (($20 == $1 || $20 == $9)) && ((index($28, $8) > 0 || index($30, $8) > 0)) && ((index($28, $16) > 0 || index($30, $16) > 0)) && ($18 == $2 || $18 == $3 || $18 == $10 || $18 == $11 || $18 + 1 == $2 || $18 + 1 == $3 || $18 + 1 == $10 || $18 + 1 == $11 || $18 - 1 == $2 || $18 - 1 == $3 || $18 - 1 == $10 || $18 - 1 == $11) && ($21 == $2 || $21 == $3 || $21 == $10 || $21 == $11 || $21 + 1 == $2 || $21 + 1 == $3 || $21 + 1 == $10 || $21 + 1 == $11 || $21 - 1 == $2 || $21 - 1 == $3 || $21 - 1 == $10 || $21 - 1 == $11)) {
		print $0
	}
}' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第四列是否匹配第一个文件的第四列
        if ($4 in file1) {
            # 如果匹配，则为第二个文件的每一行打印所有第一个文件中对应第四列的行
            for (i = 1; i <= count[$4]; i++) {
                print file1[$4][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read1_unique.bed | awk '
BEGIN {
    OFS="\t";
}
function abs(x) {
    return (x < 0) ? -x : x
}
{
		if((($39 == $1 && $39 != $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) < abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $3 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $2 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $3 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $2 + 2 && abs($40 - $21) < 1000))))) {
			if ($14 == "+") {
				$11 = $10 + 1
			} else if ($14 == "-") {
				$10 = $11 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $39, $40, $41, $42, $43, $44, $45, $9, $10, $11, $12, $13, $14, $15, $42
		} else if((($39 != $1 && $39 == $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) > abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $11 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $10 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $11 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $10 + 2 && abs($40 - $21) < 1000)))))
		{
			if ($6 == "+") {
				$3 = $2 + 1
			} else if ($6 == "-") {
				$2 = $3 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $39, $40, $41, $42, $43, $44, $45, $1, $2, $3, $4, $5, $6, $7, $42
		}
}' - > ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt
	echo "combine junction reads together"
	cat ${starres}/${fileName}Chimeric.junction.5pend.txt ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt | cut -f 15 | sort | uniq | awk -v file1=${starres}/${fileName}_read1.Chimeric.junction.5pend.txt -v file2=${starres}/${fileName}_read2.Chimeric.junction.5pend.txt -v file3=${starres}/${fileName}Chimeric.junction.5pend.txt 'BEGIN {
				FS = "\t";
        # 加载参考文件的第15列到数组
        while ((getline line < file1) > 0) {
            split(line, fields, "\t");
            part1[fields[15]] = 1;
        }
        close(file1);

        while ((getline line < file2) > 0) {
            split(line, fields, "\t");
            part2[fields[15]] = 1;
        }
        close(file2);

        while ((getline line < file3) > 0) {
            split(line, fields, "\t");
            part3[fields[15]] = 1;
        }
        close(file3);
    }
    {
        # 检查当前行是否在每个数组中
        output = $1;
        if ($1 in part1) {
            output = output "\tpart1";
        } else {
            output = output "\tnonPart1";
        }

        if ($1 in part2) {
            output = output "\tpart2";
        } else {
            output = output "\tnonPart2";
        }

        if ($1 in part3) {
            output = output "\tPE";
        } else {
            output = output "\tnonPE";
        }
        print output;
    }' > ${starres}/${fileName}_allJunctionReads
    ######
    awk 'BEGIN {FS=OFS="\t"}
NR==FNR {  # 处理 Chimeric.junction.5pend.txt 文件
    A[$15];  # 将第15列的值存储到数组A中
    next  # 跳到下一行
}
{
    # 检查条件：第二列是 nonPart1、第三列是 nonPart2、第四列是 PE
    if (($2 == "nonPart1" && $3 == "nonPart2" && $4 == "PE") ||
        ($2 == "Part1" && $3 == "Part2" && $4 == "PE")) {
        
        # 如果当前行的第2列、第3列、第4列满足条件，则检查是否第15列在数组A中
        if ($15 in A) {
            print $0  # 打印符合条件的行
        }
    }
}
' ${starres}/${fileName}Chimeric.junction.5pend.txt ${starres}/${fileName}_allJunctionReads > ${starres}/${fileName}_PE.txt
	#####
	grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$15]!=""){print}}' - ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt > ${starres}/${fileName}_read1.txt
	grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$15]!=""){print}}' - ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt > ${starres}/${fileName}_read2.txt
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($2=="nonPart1" && $3=="nonPart2" && $4=="PE"){A[$1]="yes"}else if($2=="part1" && $3=="part2" && $4=="PE"){A[$1]="yes"}}NR>FNR{if(A[$15]!=""){print}}' ${starres}/${fileName}_allJunctionReads ${starres}/${fileName}Chimeric.junction.5pend.txt > ${starres}/${fileName}_PE.txt
	grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$15]=$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$15]=$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt -| awk 'BEGIN{FS=OFS="\t"}{if($5==$11 && $8==$14){dist1=$6-$12;dist2=$9-$15;if(((dist1>=0 && dist1<=150) || (dist1<=0 && dist1>=(-150)))&&((dist2>=0 && dist2<=150) || (dist2<=0 && dist2>=(-150)))){print}}}'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$15]!=""){print}}' - ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt > ${starres}/${fileName}_read12.txt
	cat ${starres}/${fileName}_read1.txt ${starres}/${fileName}_read2.txt ${starres}/${fileName}_PE.txt ${starres}/${fileName}_read12.txt | awk 'BEGIN{FS=OFS="\t"} !seen[$1,$2,$3,$4,$8,$9,$10]++' > ${starres}/${fileName}_allJunctionReads_use_1bp
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
refstargenome=/mnt/share/Index/STAR_index/hg38_chr_149
starres=/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/1.secondbatch/8.star

for fileName in CPT21_S7_L001 ETO21_S8_L001 ETO23_S9_L001 wm-CPT-S2-R2_S3_L001 wm-ETO-S3-R2_S1_L001 wm-IR-S3-R2_S2_L001 wm-IR_S6_L001;do
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
	echo "paried-end mapping"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	rm ${starres}/${fileName}2Aligned.sortedByCoord.out.bam
	#
	echo "only consider the best chimeric and the ChimScoreDropMax<=20 and unique supported reads"
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}2Chimeric.out.junction|cut -f 10|sort|uniq -c|sed 's/^[ \t]*//g'|sed 's/ /\t/g' > ${starres}/${fileName}Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}'  ${starres}/${fileName}2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}Chim_num - > ${starres}/${fileName}Chimeric.out.junction.filter
	samtools view -h ${starres}/${fileName}Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}.chimOut.bam -tag HI -cigar > ${starres}/${fileName}.chimOut.bed
	cat ${starres}/${fileName}Chimeric.out.junction.filter | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($12 ~ /[0-9]+M/) {
        match($12, /[0-9]+M/, m);
        $21 = m[0];
    } else {
        $21 = ".";
    }
    if ($14 ~ /[0-9]+M/) {
        match($14, /[0-9]+M/, m);
        $22 = m[0];
    } else {
        $22 = ".";
    }
    print $0;
}' - > ${starres}/${fileName}Chimeric.out.junction.filter2
	cat ${starres}/${fileName}.chimOut.bed | awk '
BEGIN { OFS="\t"; }
{
    # 分割第四列，存放到第八列和第九列
    split($4, arr, "/");
    $8 = arr[1];
    $9 = arr[2];
    
    # 提取第七列中的 M 和数字
    if ($7 ~ /[0-9]+M/) {
        match($7, /[0-9]+M/, m);
        $10 = m[0];
    } else {
        $10 = ".";
    }
    
    print $0;
}' - | sort -k8,8 -k5,5n -k9,9 | awk '
BEGIN {
    OFS="\t";
}
{
    key = $8 ":" $5;
    if (key != prev_key) {
        # 如果新分组且前面有输出，则输出上一组
        if (NR > 1) {
            if (count == 2) {
                output_line = output_line "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN";
            }
            print output_line;
        }
        # 初始化新组
        output_line = $0;
        count = 1;
    } else {
        # 如果key相同，将新行串联到当前组
        count++;
        output_line = output_line "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10;
    }
    prev_key = key;
}
END {
    # 输出最后一组的数据
    if (count == 2) {
        output_line = output_line "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN" "\t" "NN";
    }
    print output_line;
}' - | awk 'BEGIN { OFS="\t"; } 
NR==FNR { aa_map[$10] = $0; next; } 
{
    if ($8 in aa_map) {
        print $0, aa_map[$8];
    } else {
        print $0, ".";
    }
}' ${starres}/${fileName}Chimeric.out.junction.filter2 - | awk '
BEGIN {
    OFS="\t";
}
function abs(x) {
    return (x < 0) ? -x : x
}
{

    # 如果第21列是 "NN"
    if ($21 == "NN") {
        if ($1 == $31 && $11 == $34 && (index($42, $10) > 0 || index($44, $10) > 0) && (index($42, $20) > 0 || index($44, $20) > 0) && (index($7, $51) > 0 || index($17, $51) > 0) && (index($7, $52) > 0 || index($17, $52) > 0) && ($32 == $2 || $32 == $12 || $32 - 1 == $3 || $32 - 1 == $13 || $32 + 1 == $3 || $32 + 1 == $13) && ($35 == $2|| $35 == $12 || $35 - 1 == $3 || $35 - 1 == $13 || $35 + 1 == $3 || $35 + 1 == $13)) {
            # 符合条件，检查第6列的符号并进行相应操作
            if ($6 == "+") {
                $3 = $2 + 1
            } else if ($6 == "-") {
                $2 = $3 - 1
            }
            if ($16 == "+") {
                $13 = $12 + 1
            } else if ($16 == "-") {
                $12 = $13 - 1
            }
            # 输出 1-7 和 11-17 列
            print $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14, $15, $16, $17, $18
        }
    }
    # 如果第21列不是 "NN"
    else {
        if ((($1 == $31 || $1 == $34) && 
             ($11 == $31 || $11 == $34) && 
             ($21 == $31 || $21 == $34) && 
             (index($42, $10) > 0 || index($44, $10) > 0) && 
             (index($42, $20) > 0 || index($44, $20) > 0) && 
             (index($42, $30) > 0 || index($44, $30) > 0) && 
             (index($7, $51) > 0 || index($17, $51) > 0 || index($27, $51) > 0) && 
             (index($7, $52) > 0 || index($17, $52) > 0 || index($27, $52) > 0) && 
             ($32 == $2 || $32 == $12 || $32 == $22 || $32 - 1 == $3 || $32 - 1 == $13 || $32 - 1 == $23 || $32 + 1 == $3 || $32 + 1 == $13 || $32 + 1 == $23) && 
             ($35 == $2 || $35 == $12 || $35 == $22 || $35 - 1 == $3 || $35 - 1 == $13 || $35 - 1 == $23 || $35 + 1 == $3 || $35 + 1 == $13 || $35 + 1 == $23))) {
            # 计算selpos
            if ($19 == 1) {
                if ($1 == $21) {
                    if ($11 == $21) {
                        if (abs($2 - $32) > abs($12 - $32)) {
                            selpos = 11
                        } else {
                            selpos = 12
                        }
                    } else {
                        selpos = 12
                    }
                } else {
                    selpos = 11
                }
            } else if ($19 == 2) {
                if ($11 == $1) {
                    if ($21 == $1) {
                        if (abs($12 - $35) > abs($22 - $35)) {
                            selpos = 22
                        } else {
                            selpos = 23
                        }
                    } else {
                        selpos = 23
                    }
                } else {
                    selpos = 22
                }
            }
            
            # 根据selpos更新列并输出不同的列
            if (selpos == 11 || selpos == 23) {
                if ($6 == "+") {
                    $3 = $2 + 1
                } else if ($6 == "-") {
                    $2 = $3 - 1
                }
                if ($26 == "+") {
                    $23 = $22 + 1
                } else if ($26 == "-") {
                    $22 = $23 - 1
                }
                # 输出1-7 和 21-27列
                print $1, $2, $3, $4, $5, $6, $7, $21, $22, $23, $24, $25, $26, $27, $28
            } else if (selpos == 12) {
                if ($16 == "+") {
                    $13 = $12 + 1
                } else if ($16 == "-") {
                    $12 = $13 - 1
                }
                if ($26 == "+") {
                    $23 = $22 + 1
                } else if ($26 == "-") {
                    $22 = $23 - 1
                }
                # 输出11-17 和 21-27列
                print $11, $12, $13, $14, $15, $16, $17, $21, $22, $23, $24, $25, $26, $27, $28
            } else if (selpos == 22) {
                if ($6 == "+") {
                    $3 = $2 + 1
                } else if ($6 == "-") {
                    $2 = $3 - 1
                }
                if ($16 == "+") {
                    $13 = $12 + 1
                } else if ($16 == "-") {
                    $12 = $13 - 1
                }
                # 输出1-7 和 11-17列
                print $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14, $15, $16, $17, $18
            }
        }
    }
}' - > ${starres}/${fileName}Chimeric.junction.5pend.txt
	echo "only consider read1"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2_read1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R1_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	echo "only consider read2"
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType WithinBAM --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1 
	STAR --runMode alignReads --genomeDir ${refstargenome} --outFileNamePrefix ${starres}/${fileName}2_read2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --outFilterMismatchNmax 999 --outFilterMultimapNmax 8 --chimSegmentMin 12 --chimScoreDropMax 400 --alignIntronMax 150 --alignMatesGapMax 150 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesIn ${trimdir}/${fileName}_R2_001.filter3.fastq --runThreadN 10 --chimOutType Junctions --alignSplicedMateMapLminOverLmate 0 --chimJunctionOverhangMin 12 --chimNonchimScoreDropMin 10 --chimMultimapNmax 8 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --chimOutJunctionFormat 1
	echo "obtain interacted read pairs for unique mapped reads"
	samtools view -h -b -q 255 ${starres}/${fileName}2_read1Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read1_unique.bam
 	samtools view -h -b -q 255 ${starres}/${fileName}2_read2Aligned.sortedByCoord.out.bam >  ${starres}/${fileName}_read2_unique.bam
 	bedtools bamtobed -i ${starres}/${fileName}_read1_unique.bam -tag HI -cigar > ${starres}/${fileName}_read1_unique.bed
 	bedtools bamtobed -i ${starres}/${fileName}_read2_unique.bam -tag HI -cigar > ${starres}/${fileName}_read2_unique.bed
 	echo "combine reads and do filteration for read1"
	awk 'BEGIN{FS=OFS="\t"} !/^#/ && $18==$19 && $16-$18<=20 {print}' ${starres}/${fileName}2_read1Chimeric.out.junction | cut -f 10 | sort | uniq -c | sed 's/^[t]*//g' | sed 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"} {n=0; for (i=1; i<=NF; i++) if ($i != "") {if (n > 0) printf "\t"; printf "%s", $i; n++} print ""}' > ${starres}/${fileName}_read1Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}' ${starres}/${fileName}2_read1Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}_read1Chim_num - > ${starres}/${fileName}_read1Chimeric.out.junction.filter
	samtools view -h ${starres}/${fileName}_read1Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}_read1.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}_read1.chimOut.bam -tag HI -cigar > ${starres}/${fileName}_read1.chimOut.bed      
	cat ${starres}/${fileName}_read1Chimeric.out.junction.filter | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($12 ~ /[0-9]+M/) {
        match($12, /[0-9]+M/, m);
        $21 = m[0];
    } else {
        $21 = ".";
    }
    if ($14 ~ /[0-9]+M/) {
        match($14, /[0-9]+M/, m);
        $22 = m[0];
    } else {
        $22 = ".";
    }
    print $0;
}' - > ${starres}/${fileName}_read1Chimeric.out.junction.filter2
	cat ${starres}/${fileName}_read1.chimOut.bed | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($7 ~ /[0-9]+M/) {
        match($7, /[0-9]+M/, m);
        $8 = m[0];
    } else {
        $8 = ".";
    }
    print $0;
}' - | sort -k4,4 -k5,5n | awk 'BEGIN {FS=OFS="\t"}
    {
        key = $4 ":" $5;
        if (key in lines) {
            lines[key] = lines[key] OFS $0;
        } else {
            lines[key] = $0;
        }
    }
    END {
        for (key in lines) {
            print lines[key];
        }
    }' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第10列是否匹配第一个文件的第四列
        if ($10 in file1) {
            for (i = 1; i <= count[$10]; i++) {
                print file1[$10][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read1Chimeric.out.junction.filter2 | awk '
BEGIN {
    OFS="\t";
}
{
	if ((($17 == $1 || $17 == $9)) && (($20 == $1 || $20 == $9)) && ((index($28, $8) > 0 || index($30, $8) > 0)) && ((index($28, $16) > 0 || index($30, $16) > 0)) && ($18 == $2 || $18 == $3 || $18 == $10 || $18 == $11 || $18 + 1 == $2 || $18 + 1 == $3 || $18 + 1 == $10 || $18 + 1 == $11 || $18 - 1 == $2 || $18 - 1 == $3 || $18 - 1 == $10 || $18 - 1 == $11) && ($21 == $2 || $21 == $3 || $21 == $10 || $21 == $11 || $21 + 1 == $2 || $21 + 1 == $3 || $21 + 1 == $10 || $21 + 1 == $11 || $21 - 1 == $2 || $21 - 1 == $3 || $21 - 1 == $10 || $21 - 1 == $11)) {
		print $0
	}
}' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第四列是否匹配第一个文件的第四列
        if ($4 in file1) {
            # 如果匹配，则为第二个文件的每一行打印所有第一个文件中对应第四列的行
            for (i = 1; i <= count[$4]; i++) {
                print file1[$4][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read2_unique.bed | awk '
BEGIN {
    OFS="\t";
}
function abs(x) {
    return (x < 0) ? -x : x
}
{
		if((($39 == $1 && $39 != $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) < abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $3 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $2 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $3 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $2 + 2 && abs($40 - $21) < 1000))))) {
			if ($14 == "+") {
				$11 = $10 + 1
			} else if ($14 == "-") {
				$10 = $11 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $9, $10, $11, $12, $13, $14, $15, $39, $40, $41, $42, $43, $44, $45, $42
		} else if((($39 != $1 && $39 == $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) > abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $11 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $10 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $11 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $10 + 2 && abs($40 - $21) < 1000)))))
			{
			if ($6 == "+") {
				$3 = $2 + 1
			} else if ($6 == "-") {
				$2 = $3 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $1, $2, $3, $4, $5, $6, $7, $39, $40, $41, $42, $43, $44, $45, $42
		}
}' - > ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt
	echo "combine reads and do filteration for read2"
	awk 'BEGIN{FS=OFS="\t"} !/^#/ && $18==$19 && $16-$18<=20 {print}' ${starres}/${fileName}2_read2Chimeric.out.junction | cut -f 10 | sort | uniq -c | sed 's/^[t]*//g' | sed 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"} {n=0; for (i=1; i<=NF; i++) if ($i != "") {if (n > 0) printf "\t"; printf "%s", $i; n++} print ""}' > ${starres}/${fileName}_read2Chim_num
	awk 'BEGIN{FS=OFS="\t"}{if($18==$19 && $16-$18<=20){print}}' ${starres}/${fileName}2_read2Chimeric.out.junction|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$2]=$1}NR>FNR{if(A[$10]==1){print $0,A[$10]}}' ${starres}/${fileName}_read2Chim_num - > ${starres}/${fileName}_read2Chimeric.out.junction.filter
	samtools view -h ${starres}/${fileName}_read2Aligned.sortedByCoord.out.bam | grep -E '(^@.*|ch:A:1)' | samtools view -Sb - > ${starres}/${fileName}_read2.chimOut.bam
	bedtools bamtobed -i ${starres}/${fileName}_read2.chimOut.bam -tag HI -cigar > ${starres}/${fileName}_read2.chimOut.bed
	cat ${starres}/${fileName}_read2Chimeric.out.junction.filter | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($12 ~ /[0-9]+M/) {
        match($12, /[0-9]+M/, m);
        $21 = m[0];
    } else {
        $21 = ".";
    }
    if ($14 ~ /[0-9]+M/) {
        match($14, /[0-9]+M/, m);
        $22 = m[0];
    } else {
        $22 = ".";
    }
    print $0;
}' - > ${starres}/${fileName}_read2Chimeric.out.junction.filter2
	cat ${starres}/${fileName}_read2.chimOut.bed | awk '
BEGIN { OFS="\t"; }
{
    # 提取第七列中的 M 和数字
    if ($7 ~ /[0-9]+M/) {
        match($7, /[0-9]+M/, m);
        $8 = m[0];
    } else {
        $8 = ".";
    }
    print $0;
}' - | sort -k4,4 -k5,5n | awk 'BEGIN {FS=OFS="\t"}
    {
        key = $4 ":" $5;
        if (key in lines) {
            lines[key] = lines[key] OFS $0;
        } else {
            lines[key] = $0;
        }
    }
    END {
        for (key in lines) {
            print lines[key];
        }
    }' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第10列是否匹配第一个文件的第四列
        if ($10 in file1) {
            for (i = 1; i <= count[$10]; i++) {
                print file1[$10][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read2Chimeric.out.junction.filter2 | awk '
BEGIN {
    OFS="\t";
}
{
	if ((($17 == $1 || $17 == $9)) && (($20 == $1 || $20 == $9)) && ((index($28, $8) > 0 || index($30, $8) > 0)) && ((index($28, $16) > 0 || index($30, $16) > 0)) && ($18 == $2 || $18 == $3 || $18 == $10 || $18 == $11 || $18 + 1 == $2 || $18 + 1 == $3 || $18 + 1 == $10 || $18 + 1 == $11 || $18 - 1 == $2 || $18 - 1 == $3 || $18 - 1 == $10 || $18 - 1 == $11) && ($21 == $2 || $21 == $3 || $21 == $10 || $21 == $11 || $21 + 1 == $2 || $21 + 1 == $3 || $21 + 1 == $10 || $21 + 1 == $11 || $21 - 1 == $2 || $21 - 1 == $3 || $21 - 1 == $10 || $21 - 1 == $11)) {
		print $0
	}
}' - | awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {
        # 读取第一个文件，将第四列作为键，每行内容存为数组的值
        file1[$4][++count[$4]] = $0;
        next;
    }
    {
        # 遍历第二个文件，检查第四列是否匹配第一个文件的第四列
        if ($4 in file1) {
            # 如果匹配，则为第二个文件的每一行打印所有第一个文件中对应第四列的行
            for (i = 1; i <= count[$4]; i++) {
                print file1[$4][i], $0;
            }
        }
    }' - ${starres}/${fileName}_read1_unique.bed | awk '
BEGIN {
    OFS="\t";
}
function abs(x) {
    return (x < 0) ? -x : x
}
{
		if((($39 == $1 && $39 != $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) < abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $3 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $2 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $3 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $2 + 2 && abs($40 - $21) < 1000))))) {
			if ($14 == "+") {
				$11 = $10 + 1
			} else if ($14 == "-") {
				$10 = $11 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $39, $40, $41, $42, $43, $44, $45, $9, $10, $11, $12, $13, $14, $15, $42
		} else if((($39 != $1 && $39 == $9) || ($39 == $1 && $39 == $9 && abs($2 - $40) > abs($10 - $40))) && (($39 == $17 && (($18 > $41 - 2 && $18 > $11 - 2 && abs($18 - $41) < 1000) || ($18 < $40 + 2 && $18 < $10 + 2 && abs($40 - $18) < 1000))) || ($39 == $20 && (($21 > $41 - 2 && $21 > $11 - 2 && abs($21 - $41) < 1000) || ($21 < $40 + 2 && $21 < $10 + 2 && abs($40 - $21) < 1000)))))
		{
			if ($6 == "+") {
				$3 = $2 + 1
			} else if ($6 == "-") {
				$2 = $3 - 1
			}
			if ($44 == "+") {
				$41 = $40 + 1
			} else if ($44 == "-") {
				$40 = $41 - 1
			}
			# 输出9-15 和 39-45列
			print $39, $40, $41, $42, $43, $44, $45, $1, $2, $3, $4, $5, $6, $7, $42
		}
}' - > ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt
	echo "combine junction reads together"
	cat ${starres}/${fileName}Chimeric.junction.5pend.txt ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt | cut -f 15 | sort | uniq | awk -v file1=${starres}/${fileName}_read1.Chimeric.junction.5pend.txt -v file2=${starres}/${fileName}_read2.Chimeric.junction.5pend.txt -v file3=${starres}/${fileName}Chimeric.junction.5pend.txt 'BEGIN {
				FS = "\t";
        # 加载参考文件的第15列到数组
        while ((getline line < file1) > 0) {
            split(line, fields, "\t");
            part1[fields[15]] = 1;
        }
        close(file1);

        while ((getline line < file2) > 0) {
            split(line, fields, "\t");
            part2[fields[15]] = 1;
        }
        close(file2);

        while ((getline line < file3) > 0) {
            split(line, fields, "\t");
            part3[fields[15]] = 1;
        }
        close(file3);
    }
    {
        # 检查当前行是否在每个数组中
        output = $1;
        if ($1 in part1) {
            output = output "\tpart1";
        } else {
            output = output "\tnonPart1";
        }

        if ($1 in part2) {
            output = output "\tpart2";
        } else {
            output = output "\tnonPart2";
        }

        if ($1 in part3) {
            output = output "\tPE";
        } else {
            output = output "\tnonPE";
        }
        print output;
    }' > ${starres}/${fileName}_allJunctionReads
    ######
    awk 'BEGIN {FS=OFS="\t"}
NR==FNR {  # 处理 Chimeric.junction.5pend.txt 文件
    A[$15];  # 将第15列的值存储到数组A中
    next  # 跳到下一行
}
{
    # 检查条件：第二列是 nonPart1、第三列是 nonPart2、第四列是 PE
    if (($2 == "nonPart1" && $3 == "nonPart2" && $4 == "PE") ||
        ($2 == "Part1" && $3 == "Part2" && $4 == "PE")) {
        
        # 如果当前行的第2列、第3列、第4列满足条件，则检查是否第15列在数组A中
        if ($15 in A) {
            print $0  # 打印符合条件的行
        }
    }
}
' ${starres}/${fileName}Chimeric.junction.5pend.txt ${starres}/${fileName}_allJunctionReads > ${starres}/${fileName}_PE.txt
	#####
	grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'nonPart2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$15]!=""){print}}' - ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt > ${starres}/${fileName}_read1.txt
	grep -w 'nonPart1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$15]!=""){print}}' - ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt > ${starres}/${fileName}_read2.txt
	awk 'BEGIN{FS=OFS="\t"}NR==FNR{if($2=="nonPart1" && $3=="nonPart2" && $4=="PE"){A[$1]="yes"}else if($2=="part1" && $3=="part2" && $4=="PE"){A[$1]="yes"}}NR>FNR{if(A[$15]!=""){print}}' ${starres}/${fileName}_allJunctionReads ${starres}/${fileName}Chimeric.junction.5pend.txt > ${starres}/${fileName}_PE.txt
	grep -w 'part1' ${starres}/${fileName}_allJunctionReads |grep -w 'part2'|grep -w 'nonPE'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$15]=$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt -|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$15]=$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10}NR>FNR{print $0,A[$1]}' ${starres}/${fileName}_read2.Chimeric.junction.5pend.txt -| awk 'BEGIN{FS=OFS="\t"}{if($5==$11 && $8==$14){dist1=$6-$12;dist2=$9-$15;if(((dist1>=0 && dist1<=150) || (dist1<=0 && dist1>=(-150)))&&((dist2>=0 && dist2<=150) || (dist2<=0 && dist2>=(-150)))){print}}}'|awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]="yes"}NR>FNR{if(A[$15]!=""){print}}' - ${starres}/${fileName}_read1.Chimeric.junction.5pend.txt > ${starres}/${fileName}_read12.txt
	cat ${starres}/${fileName}_read1.txt ${starres}/${fileName}_read2.txt ${starres}/${fileName}_PE.txt ${starres}/${fileName}_read12.txt | awk 'BEGIN{FS=OFS="\t"} !seen[$1,$2,$3,$4,$8,$9,$10]++' > ${starres}/${fileName}_allJunctionReads_use_1bp
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



########################################
######################################## combine two replicates for MCF-7 without drug treatment
########################################

R CMD BATCH ./2.combine1.r
## outputs: merged bait sites, translocation sites after removing duplications


