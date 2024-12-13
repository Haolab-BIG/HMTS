#! /bin/bash

################################################ firstbatch
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
        if ($1 == $31 && $11 == $34 && (index($42, $10) > 0 || index($44, $10) > 0) && (index($42, $20) > 0 || index($44, $20) > 0) && (index($7, $51) > 0 || index($17, $51) > 0) && (index($7, $52) > 0 || index($17, $52) > 0)) {
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
             (index($7, $52) > 0 || index($17, $52) > 0 || index($27, $52) > 0))) {
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
done






