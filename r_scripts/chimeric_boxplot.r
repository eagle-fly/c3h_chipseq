
threshold <- "stringent"
#threshold <- "relaxed"

setwd(paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, sep=""))
infile <- "sorted_chim/chim_all.txt"
chim_data <- read.table(infile, header=TRUE, sep="\t")
load("/Users/ying/Work/BCCRC/Projects/colon/artem_scripts/RM/RepeatMasker.Rdata")


TE_total <- 1407573768
LINE_total <- claStats[3,3]
SINE_total <- claStats[15,3]
LTR_total <- claStats[6,3]

TE_ratio <- TE_total / TE_total
LINE_ratio <- LINE_total / TE_total
SINE_ratio <- SINE_total / TE_total
LTR_ratio <- LTR_total / TE_total

TE_normal <- chim_data[chim_data$LIBRARY%%2==1,]
TE_nor_bySample <- as.vector(table(TE_normal$LIBRARY))
TE_chim_mean <- mean(TE_nor_bySample)
TE_cancer <- chim_data[chim_data$LIBRARY%%2==0,]
TE_can_bySample <- as.vector(table(TE_cancer$LIBRARY))

LINE_normal <- chim_data[chim_data$repClass=="LINE" & chim_data$LIBRARY%%2==1,]
LINE_cancer <- chim_data[chim_data$repClass=="LINE" & chim_data$LIBRARY%%2==0,]
SINE_normal <- chim_data[chim_data$repClass=="SINE" & chim_data$LIBRARY%%2==1,]
SINE_cancer <- chim_data[chim_data$repClass=="SINE" & chim_data$LIBRARY%%2==0,]
LTR_normal <- chim_data[chim_data$repClass=="LTR" & chim_data$LIBRARY%%2==1,]
LTR_cancer <- chim_data[chim_data$repClass=="LTR" & chim_data$LIBRARY%%2==0,]

LINE_nor_bySample <- table(LINE_normal$LIBRARY)
LINE_can_bySample <- table(LINE_cancer$LIBRARY)
SINE_nor_bySample <- table(SINE_normal$LIBRARY)
SINE_can_bySample <- table(SINE_cancer$LIBRARY)
LTR_nor_bySample <- table(LTR_normal$LIBRARY)
LTR_can_bySample <- table(LTR_cancer$LIBRARY)

LINE_nor_corrected <- as.vector(LINE_nor_bySample) / TE_chim_mean / LINE_ratio 
LINE_can_corrected <- as.vector(LINE_can_bySample) / TE_chim_mean / LINE_ratio
SINE_nor_corrected <- as.vector(SINE_nor_bySample) / TE_chim_mean / SINE_ratio
SINE_can_corrected <- as.vector(SINE_can_bySample) / TE_chim_mean / SINE_ratio
LTR_nor_corrected <- as.vector(LTR_nor_bySample) / TE_chim_mean / LTR_ratio
LTR_can_corrected <- as.vector(LTR_can_bySample) / TE_chim_mean / LTR_ratio


# boxplot of overall chimera expression of all TEs pooled
y_lim <- max(c(TE_nor_bySample, TE_can_bySample))
png(paste("../../results/", threshold, "/chimExp_pooled.png", sep=""), width=6, height=6, units="in", res=300)
par(mgp=c(4,1,0), mar=c(6,6,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("orange", "blue")
boxplot(list(TE_nor_bySample, TE_can_bySample), range=1, col=colorset,
        ylab="Chimeric Transcripts in Sample", 
        cex.axis=0.9, cex.lab=1.2, names=c("Normal", "Cancer"))
dev.off()


# boxplot of overall chimera expression of major TE classes
y_lim <- max(c(LINE_nor_corrected, LINE_can_corrected, SINE_nor_corrected, SINE_can_corrected, LTR_nor_corrected, LTR_can_corrected))
png(paste("../../results/", threshold, "/chimExp_all.png", sep=""), width=6, height=6, units="in", res=300)
par(mgp=c(4,1,0), mar=c(6,6,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("orange", "blue")
boxplot(list(LINE_nor_corrected, LINE_can_corrected, SINE_nor_corrected, SINE_can_corrected, 
             LTR_nor_corrected, LTR_can_corrected), range=1, col=colorset,
            main="All Chimeric", xlab="Repeat Class", ylab="Normalized Chimeric Transcripts", 
            names=c("        LINE", "", "        SINE", "", "        LTR", ""), cex.axis=0.9, cex.lab=1.2,
            at=c(1,2,4,5,7,8))
legend(1, y_lim, c("Normal", "Cancer"), colorset)
segments(0,1,9,1, col="red", lty=2)
dev.off()


#================================================================================

# barplot of chimera expression comparison between normal and cancer for each TE class
listfile <- "../list.txt"
sample_list <- read.table (listfile, header=FALSE)[,1]
LINE_can_index <- which(sample_list %in% rownames(LINE_can_bySample))
LINE_nor_index <- which(sample_list %in% rownames(LINE_nor_bySample))
LINE_bySample <- rep(0, length(sample_list))
LINE_bySample[LINE_can_index] <- LINE_can_bySample
LINE_bySample[LINE_nor_index] <- LINE_nor_bySample

png(paste("../../results/", threshold, "/chimExp_bar_LINE.png", sep=""), width=16, height=8, units="in", res=300)
par(mgp=c(3,1,0), mar=c(5,5,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("darkblue", "orange")
barplot(matrix(LINE_bySample, 2, length(sample_list)/2), beside=TRUE, col=colorset, xlab="Samples", ylab="Number of Chimeras")
axis(1, tick=FALSE, labels=FALSE)
text(x=sort(c(seq(1,187,3), seq(2,188,3)))+0.3, y=par("usr")[3]-0.2, srt=90, adj = 0.8, labels = sample_list, cex=0.5, xpd = TRUE)
y_lim <- max(c(LINE_bySample))
legend(1, y_lim, c("Cancer", "Normal"), colorset)
dev.off()


SINE_can_index <- which(sample_list %in% rownames(SINE_can_bySample))
SINE_nor_index <- which(sample_list %in% rownames(SINE_nor_bySample))
SINE_bySample <- rep(0, length(sample_list))
SINE_bySample[SINE_can_index] <- SINE_can_bySample
SINE_bySample[SINE_nor_index] <- SINE_nor_bySample

png(paste("../../results/", threshold, "/chimExp_bar_SINE.png", sep=""), width=16, height=8, units="in", res=300)
par(mgp=c(3,1,0), mar=c(5,5,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("darkblue", "orange")
barplot(matrix(SINE_bySample, 2, length(sample_list)/2), beside=TRUE, col=colorset, xlab="Samples", ylab="Number of Chimeras")
axis(1, tick=FALSE, labels=FALSE)
text(x=sort(c(seq(1,187,3), seq(2,188,3)))+0.3, y=par("usr")[3]-0.2, srt=90, adj = 0.8, labels = sample_list, cex=0.5, xpd = TRUE)
y_lim <- max(c(SINE_bySample))
legend(1, y_lim, c("Cancer", "Normal"), colorset)
dev.off()


LTR_can_index <- which(sample_list %in% rownames(LTR_can_bySample))
LTR_nor_index <- which(sample_list %in% rownames(LTR_nor_bySample))
LTR_bySample <- rep(0, length(sample_list))
LTR_bySample[LTR_can_index] <- LTR_can_bySample
LTR_bySample[LTR_nor_index] <- LTR_nor_bySample

png(paste("../../results/", threshold, "/chimExp_bar_LTR.png", sep=""), width=16, height=8, units="in", res=300)
par(mgp=c(3,1,0), mar=c(5,5,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("darkblue", "orange")
barplot(matrix(LTR_bySample, 2, length(sample_list)/2), beside=TRUE, col=colorset, xlab="Samples", ylab="Number of Chimeras")
axis(1, tick=FALSE, labels=FALSE)
text(x=sort(c(seq(1,187,3), seq(2,188,3)))+0.3, y=par("usr")[3]-0.2, srt=90, adj = 0.8, labels = sample_list, cex=0.5, xpd = TRUE)
y_lim <- max(c(LTR_bySample))
legend(1, y_lim, c("Cancer", "Normal"), colorset)
dev.off()


#================================================================================

repFamilies <- c("L1", "L2", "Other")
repFamSize <- c(famStats[24,4], famStats[26,4], famStats[4,4]+famStats[8,4]+famStats[39,4]+famStats[40,4])
total_families <- length(repFamilies)
chimValues <- list();
for (i in 1:total_families) {
  normal <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==1,]
  cancer <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==0,]
  repFam_ratio <- repFamSize[i] / TE_total
  nor_bySample <- table(normal$LIBRARY)
  nor_bySample[nor_bySample < 3] <- 0 # if the # of cases is less than 3, don't count
  null_samples <- 63 - length(nor_bySample)
  nor_bySample <- c(nor_bySample, rep(0, null_samples))
  nor_bySample <- as.vector(nor_bySample) / TE_chim_mean / repFam_ratio
  can_bySample <- table(cancer$LIBRARY)
  can_bySample[can_bySample < 3] <- 0 # if the # of cases is less than 3, don't count
  null_samples <- 63 - length(can_bySample)
  can_bySample <- c(can_bySample, rep(0, null_samples))
  can_bySample <- as.vector(can_bySample) / TE_chim_mean / repFam_ratio
  chimValues <- append(chimValues, list(nor_bySample, can_bySample))
}

y_lim <- max(unlist(chimValues))
png(paste("../../results/", threshold, "/chimExp_LINE.png", sep=""), width=6, height=6, units="in", res=300)
colorset <- c("orange", "blue")
boxplot(chimValues, range=1, col=colorset,
        main="LINE Chimeric", xlab="Repeat Family", ylab="Normalized Chimeric Transcripts", 
        names=c("        L1", "", "        L2", "", "        Others", ""), cex.axis=0.9, cex.lab=1.2,
        at=c(1,2,4,5,7,8))
legend(0.5, y_lim, c("Normal", "Cancer"), colorset)
segments(0,1,9,1, col="red", lty=2)
dev.off()

#====
# SINEs
repFamilies <- c("Alu", "MIR", "Other")
repFamSize <- c(famStats[2,4], famStats[31,4], famStats[5,4]+famStats[44,4])
total_families <- length(repFamilies)
chimValues <- list();
for (i in 1:total_families) {
  normal <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==1,]
  cancer <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==0,]
  repFam_ratio <- repFamSize[i] / TE_total
  nor_bySample <- table(normal$LIBRARY)
  nor_bySample[nor_bySample < 3] <- 0 # if the # of cases is less than 3, don't count
  null_samples <- 63 - length(nor_bySample)
  nor_bySample <- c(nor_bySample, rep(0, null_samples))
  nor_bySample <- as.vector(nor_bySample) / TE_chim_mean / repFam_ratio
  can_bySample <- table(cancer$LIBRARY)
  can_bySample[can_bySample < 3] <- 0 # if the # of cases is less than 3, don't count
  null_samples <- 63 - length(can_bySample)
  can_bySample <- c(can_bySample, rep(0, null_samples))
  can_bySample <- as.vector(can_bySample) / TE_chim_mean / repFam_ratio
  chimValues <- append(chimValues, list(nor_bySample, can_bySample))
}

y_lim <- max(unlist(chimValues))
png(paste("../../results/", threshold, "/chimExp_SINE.png", sep=""), width=6, height=6, units="in", res=300)
colorset <- c("orange", "blue")
boxplot(chimValues, range=1, col=colorset,
        main="SINE Chimeric", xlab="Repeat Family", ylab="Normalized Chimeric Transcripts", 
        names=c("        Alu", "", "        MIR", "", "        Others", ""), cex.axis=0.9, cex.lab=1.2,
        at=c(1,2,4,5,7,8))
legend(0.5, y_lim, c("Normal", "Cancer"), colorset)
segments(0,1,9,1, col="red", lty=2)
dev.off()

#====
# LTRs

repFamilies <- c("ERV1", "ERVL", "ERVL-MaLR", "Other")
repFamSize <- c(famStats[10,4], famStats[12,4], famStats[14,4], 
                famStats[9,4]+famStats[11,4]+famStats[13,4]+famStats[15,4]+famStats[16,4]+famStats[28,4])
total_families <- length(repFamilies)
chimValues <- list();
for (i in 1:total_families) {
  normal <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==1,]
  cancer <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==0,]
  repFam_ratio <- repFamSize[i] / TE_total
  nor_bySample <- table(normal$LIBRARY)
  nor_bySample[nor_bySample < 3] <- 0 # if the # of cases is less than 3, don't count
  null_samples <- 63 - length(nor_bySample)
  nor_bySample <- c(nor_bySample, rep(0, null_samples))
  nor_bySample <- as.vector(nor_bySample) / TE_chim_mean / repFam_ratio
  can_bySample <- table(cancer$LIBRARY)
  can_bySample[can_bySample < 3] <- 0 # if the # of cases is less than 3, don't count
  null_samples <- 63 - length(can_bySample)
  can_bySample <- c(can_bySample, rep(0, null_samples))
  can_bySample <- as.vector(can_bySample) / TE_chim_mean / repFam_ratio
  chimValues <- append(chimValues, list(nor_bySample, can_bySample))
}

y_lim <- max(unlist(chimValues))
png(paste("../../results/", threshold, "/chimExp_LTR.png", sep=""), width=6, height=6, units="in", res=300)
colorset <- c("orange", "blue")
boxplot(chimValues, range=1, col=colorset,
        main="LTR Chimeric", xlab="Repeat Family", ylab="Normalized Chimeric Transcripts", 
        names=c("        ERV1", "", "        ERVL", "", "        MaLR", "", "        Others", ""), cex.axis=0.9, cex.lab=1.2,
        at=c(1,2,4,5,7,8,10,11))
legend(9, y_lim, c("Normal", "Cancer"), colorset)
segments(0,1,12,1, col="red", lty=2)
dev.off()
