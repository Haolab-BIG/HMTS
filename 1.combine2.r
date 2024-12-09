##########
options(stringsAsFactors=FALSE)
options(scipen = 999)
setwd("/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/4.junctionreads")
###
baitAndTransfiles <- list.files(path="/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/3.star", pattern="baitAndTrans", full.names=T)
DOX1 <- read.table(baitAndTransfiles[1], header=F, sep="\t", quote="")
DOX2 <- read.table(baitAndTransfiles[2], header=F, sep="\t", quote="")
DOXcombine <- rbind(DOX1, DOX2)

##########
GetR1R2bed <- function(combineJR, name, distance=1000){
	combineJR <- combineJR[combineJR$V1 %in% paste0("chr", c(1:22, "X", "Y")) & combineJR$V4 %in% paste0("chr", c(1:22, "X", "Y")), ]
	BaitR1 <- combineJR[, c(1, 2, 2, 7, 7, 3)]
	BaitR1[, 2] <- BaitR1[, 2]-1
	BaitR1[, 5] <- "."
	BaitR2 <- combineJR[, c(4, 5, 5, 7, 7, 6)]
	BaitR2[, 2] <- BaitR2[, 2]-1
	BaitR2[, 5] <- "."
	write.table(BaitR1, file=paste0("Bait_", name, ".bed"), col.names=F, row.names=F, sep="\t", quote=F)
	write.table(BaitR2, file=paste0("Translocation_", name, ".bed"), col.names=F, row.names=F, sep="\t", quote=F)
	system(paste0("sort -V -k1,1 -k2,2 ", "Bait_", name, ".bed", " > ", "Bait_", name, ".sorted.bed"))
	#system(paste0("sort -V -k1,1 -k2,2 ", "Translocation_", name, ".bed", " > ", "Translocation_", name, ".sorted.bed"))
	## merge bait
	system(paste0("bedtools merge -i ", "Bait_", name, ".sorted.bed", " -d ", distance, " -c 4,5,6 -o collapse,distinct,distinct", " > ", "Bait_", name, ".merged1kb.bed"))
	baitmerge <- read.table(paste0("Bait_", name, ".merged1kb.bed"), header=F, sep="\t", quote="")
	baitmerge$V4 <- paste0(name, ".merged1kb", "_bait_", 1:(dim(baitmerge)[1]))
	write.table(baitmerge[, 1:4], file=paste0("Bait_", name, ".merged1kb2.bed"), col.names=F, row.names=F, sep="\t", quote=F)
	print(dim(baitmerge))
	system(paste0("sort -V -k1,1 -k2,2 ", "Bait_", name, ".merged1kb2.bed", " > ", "Bait_", name, ".merged1kb2_sort.bed"))
	system(paste0("bedtools intersect -wa -wb -a ", "Bait_", name, ".merged1kb2_sort.bed", " -b ", "Bait_", name, ".sorted.bed", " -sorted > ", "Bait_", name, ".mapp.bed"))
	## name translocation
	translocationmerge <- BaitR2
	colnames(translocationmerge) <- paste0("V", 1:length(colnames(translocationmerge)))
	translocationmerge$V7 <- paste0(translocationmerge[, 1], ":", translocationmerge[, 2], ":", translocationmerge[, 3])
	translocationdf <- data.frame(loc=unique(translocationmerge$V7), name=paste0(name, "_translocation_", 1:(length(unique(translocationmerge$V7)))))
	translocationdf2 <- merge(translocationmerge, translocationdf, by.x="V7", by.y="loc")
	translocationmerge2 <- translocationdf2[, c(2:4, 8, 6:7, 2:7)]
	translocationdf3 <- data.frame(do.call(rbind, strsplit(translocationdf[, 1], "\\:")), translocationdf[, 2])  ## 4col
	print(dim(translocationdf3)[1])
	write.table(translocationdf3, file=paste0("Translocation_", name, "2.bed"), col.names=F, row.names=F, sep="\t", quote=F)
	system(paste0("sort -V -k1,1 -k2,2 ", "Translocation_", name, "2.bed", " > ", "Translocation_", name, "2_sort.bed"))
	## find relationship between bait and translocation
	bb <- read.table(paste0("Bait_", name, ".mapp.bed"), header=F, sep="\t", quote="")[, c(4, 8)]
	tt <- translocationmerge2[, c(4, 10)]
	mm <- merge(bb, tt, by.x="V8", by.y="V4")
	uu <- unique(mm[, 2:3])
	print(dim(uu)[1])
	write.table(uu, file=paste0(name, "_btunique2.txt"), col.names=F, row.names=F, sep="\t", quote=F)
}

###
GetR1R2bed(combineJR=DOX1, name="DOX1")
GetR1R2bed(combineJR=DOX2, name="DOX2")
GetR1R2bed(combineJR=DOXcombine, name="DOXcombine")
