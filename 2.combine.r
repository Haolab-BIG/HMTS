####### step1: combine all reads supporting baits and translocations
options(stringsAsFactors=FALSE)
options(scipen = 999)
setwd("/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/11.junctionreads")
#
print("step1: combine all reads supporting baits and translocations")
alllocfile <- list.files(path="/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/10.star", pattern="_allBaitTransReads_1bp", full.names=T)

####### step2: find pairs and merge bait within 1kb, remove duplications
print("step2: find pairs and merge bait within 1kb, remove duplications")

GetR1R2bed <- function(combineJR, name, distance1=1000, distance2=1000){
	BaitR1 <- combineJR[, c(1:3, 9, 9, 4)]
	BaitR1[, 5] <- 1
	BaitR2 <- combineJR[, c(5:7, 9, 9, 8)]
	BaitR2[, 5] <- 1
	write.table(BaitR1, file=paste0("Bait_", name, ".bed"), col.names=F, row.names=F, sep="\t", quote=F)
	write.table(BaitR2, file=paste0("Translocation_", name, ".bed"), col.names=F, row.names=F, sep="\t", quote=F)
	system(paste0("sort -V -k1,1 -k2,2 ", "Bait_", name, ".bed", " > ", "Bait_", name, ".sorted.bed"))
	system(paste0("sort -V -k1,1 -k2,2 ", "Translocation_", name, ".bed", " > ", "Translocation_", name, ".sorted.bed"))
	
	## merge bait
	
	print("merge bait")
	if(distance1==1000 & distance2==1000){
		system(paste0("bedtools merge -i ", "Bait_", name, ".sorted.bed", " -d ", distance1, " > ", "Bait_", name, ".merged1kb.bed"))
		system(paste0("bedtools merge -i ", "Translocation_", name, ".sorted.bed", " -d ", distance2, " > ", "Translocation_", name, ".merged1kb.bed"))
	} else {
		system(paste0("bedtools merge -i ", "Bait_", name, ".sorted.bed", " -d ", distance, " > ", "Bait_", name, ".merged.", distance1, ".bed"))
		system(paste0("bedtools merge -i ", "Translocation_", name, ".sorted.bed", " -d ", distance, " > ", "Translocation_", name, ".merged.", distance2, ".bed"))
	}
	baitmerge <- read.table(paste0("Bait_", name, ".merged1kb.bed"), header=F, sep="\t", quote="")
	baitmerge$V4 <- paste0(name, ".merged1kb", "_bait_", 1:(dim(baitmerge)[1]))
	write.table(baitmerge[, 1:4], file=paste0("Bait_", name, ".merged1kb2.bed"), col.names=F, row.names=F, sep="\t", quote=F)
	print(dim(baitmerge))
	system(paste0("sort -V -k1,1 -k2,2 ", "Bait_", name, ".merged1kb2.bed", " > ", "Bait_", name, ".merged1kb2_sort.bed"))
	system(paste0("bedtools intersect -wa -wb -a ", "Bait_", name, ".merged1kb2_sort.bed", " -b ", "Bait_", name, ".sorted.bed", " -sorted > ", "Bait_", name, ".mapp.bed"))

	## no merge translocation
	
	print("no merge translocation")
	translocationmerge <- read.table(paste0("Translocation_", name, ".sorted.bed"), header=F, sep="\t", quote="")
	translocationmerge$V7 <- paste0(translocationmerge[, 1], ":", translocationmerge[, 2], ":", translocationmerge[, 3])
	translocationdf <- data.frame(loc=unique(translocationmerge$V7), name=paste0(name, "_translocation_", 1:(length(unique(translocationmerge$V7)))))
	translocationdf2 <- merge(translocationmerge, translocationdf, by.x="V7", by.y="loc")
	translocationmerge2 <- translocationdf2[, c(2:4, 8, 6:7, 2:7)]
	translocationdf3 <- data.frame(do.call(rbind, strsplit(translocationdf[, 1], "\\:")), translocationdf[, 2])  ## 4col
	print(dim(translocationdf3)[1])
	write.table(translocationdf3, file=paste0("Translocation_", name, "2.bed"), col.names=F, row.names=F, sep="\t", quote=F)
	system(paste0("sort -V -k1,1 -k2,2 ", "Translocation_", name, "2.bed", " > ", "Translocation_", name, "2_sort.bed"))	
	bb <- read.table(paste0("Bait_", name, ".mapp.bed"), header=F, sep="\t", quote="")[, c(4, 8)]
	tt <- translocationmerge2[, c(4, 10)]
	mm <- merge(bb, tt, by.x="V8", by.y="V4")
	uu <- unique(mm[, 2:3])
	print(dim(uu)[1])
	write.table(uu, file=paste0(name, "_btunique2.txt"), col.names=F, row.names=F, sep="\t", quote=F)
	
	## merge translocation
	
	print("merge translocation")
	transmerge <- read.table(paste0("Translocation_", name, ".merged1kb.bed"), header=F, sep="\t", quote="")
	transmerge$V4 <- paste0(name, ".merged1kb", "Translocation_", 1:(dim(transmerge)[1]))
	write.table(transmerge[, 1:4], file=paste0("Translocation_", name, ".merged1kb2.bed"), col.names=F, row.names=F, sep="\t", quote=F)
	#print(dim(transmerge))
	system(paste0("sort -V -k1,1 -k2,2 ", "Translocation_", name, ".merged1kb2.bed", " > ", "Translocation_", name, ".merged1kb2_sort.bed"))
	system(paste0("bedtools intersect -wa -wb -a ", "Translocation_", name, ".merged1kb2_sort.bed", " -b ", "Translocation_", name, ".sorted.bed", " -sorted > ", "Translocation_", name, ".mapp.bed"))
	bb2 <- read.table(paste0("Bait_", name, ".mapp.bed"), header=F, sep="\t", quote="")[, c(4, 8)]
	tt2 <- read.table(paste0("Translocation_", name, ".mapp.bed"), header=F, sep="\t", quote="")[, c(4, 8)]
	mm2 <- merge(bb2, tt2, by="V8")
	uu2 <- unique(mm2[, 2:3])
	#print(dim(uu2)[1])
	write.table(uu2, file=paste0(name, "_btunique.txt"), col.names=F, row.names=F, sep="\t", quote=F)

	## no merge bait
	print("no merge bait")
	BaitR1$V7 <- paste0(BaitR1[, 1], ":", BaitR1[, 2], ":", BaitR1[, 3])
	baitdf <- data.frame(loc=unique(BaitR1$V7), name=paste0(name, "_Bait_", 1:(length(unique(BaitR1$V7)))))
	baitdf2 <- merge(BaitR1, baitdf, by.x="V7", by.y="loc")
	baitmerge2 <- baitdf2[, c(2:4, 8, 6:7, 2:7)]
	baitdf3 <- data.frame(do.call(rbind, strsplit(baitdf[, 1], "\\:")), baitdf[, 2])  ## 4col
	print(dim(baitdf3)[1])
	write.table(baitdf3, file=paste0("Bait_", name, "2.bed"), col.names=F, row.names=F, sep="\t", quote=F)
	system(paste0("sort -V -k1,1 -k2,2 ", "Bait_", name, "2.bed", " > ", "Bait_", name, "2_sort.bed"))	
	bb3 <- baitmerge2[, c(4, 10)]
	tt3 <- translocationmerge2[, c(4, 10)]
	mm3 <- merge(bb3, tt3, by.x="V9", by.y="V4")
	uu3 <- unique(mm3[, 2:3])
	print(dim(uu3)[1])
	write.table(uu3, file=paste0(name, "_btunique3.txt"), col.names=F, row.names=F, sep="\t", quote=F)
}

for(ii in 1:length(alllocfile)){
	print(alllocfile[ii])
	namefort <- gsub("/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/10.star/", "", gsub("_allBaitTransReads_1bp", "", alllocfile[ii]))
	JRinfor <- read.table(alllocfile[ii], header=F, sep="\t", quote="")
	GetR1R2bed(combineJR=JRinfor, name=namefort)
}


aa <- read.table("/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/10.star/wm-DOX1_S4_L001_allBaitTransReads_1bp", header=F, sep="\t", quote="")
bb <- read.table("/mnt1/2.NAS2024/wutan/4.translocation/Degron/5.DOX/10.star/wm-DOX2_S5_L001_allBaitTransReads_1bp", header=F, sep="\t", quote="")
JRinfor <- rbind(aa, bb)
namefort <- "DOXcombine"
print(namefort)
GetR1R2bed(combineJR=JRinfor, name=namefort)
