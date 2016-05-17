infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/sorted_chim_relaxed/chim_all.txt"
chim_data <- read.table(infile, header=TRUE, sep="\t")


LINE_total <- 1570523
SINE_total <- 1852545
LTR_total <- 748597

LINE_ratio <- LINE_total / LTR_total
SINE_ratio <- SINE_total / LTR_total
LTR_ratio <- LTR_total / LTR_total


LINE_normal <- chim_data[chim_data$repClass=="LINE" & chim_data$LIBRARY%%2==1,]
LINE_cancer <- chim_data[chim_data$repClass=="LINE" & chim_data$LIBRARY%%2==0,]
SINE_normal <- chim_data[chim_data$repClass=="SINE" & chim_data$LIBRARY%%2==1,]
SINE_cancer <- chim_data[chim_data$repClass=="SINE" & chim_data$LIBRARY%%2==0,]
LTR_normal <- chim_data[chim_data$repClass=="LTR" & chim_data$LIBRARY%%2==1,]
LTR_cancer <- chim_data[chim_data$repClass=="LTR" & chim_data$LIBRARY%%2==0,]

LINE_nor_bySample <- table(LINE_normal$LIBRARY) / LINE_ratio
LINE_can_bySample <- table(LINE_cancer$LIBRARY) / LINE_ratio
SINE_nor_bySample <- table(SINE_normal$LIBRARY) / SINE_ratio
SINE_can_bySample <- table(SINE_cancer$LIBRARY) / SINE_ratio
LTR_nor_bySample <- table(LTR_normal$LIBRARY) / LTR_ratio
LTR_can_bySample <- table(LTR_cancer$LIBRARY) / LTR_ratio

LINE_nor_mean <- mean(LINE_nor_bySample)
LINE_can_mean <- mean(LINE_can_bySample)
SINE_nor_mean <- mean(SINE_nor_bySample)
SINE_can_mean <- mean(SINE_can_bySample)
LTR_nor_mean <- mean(LTR_nor_bySample)
LTR_can_mean <- mean(LTR_can_bySample)
LINE_nor_sd <- sd(LINE_nor_bySample)
LINE_can_sd <- sd(LINE_can_bySample)
SINE_nor_sd <- sd(SINE_nor_bySample)
SINE_can_sd <- sd(SINE_can_bySample)
LTR_nor_sd <- sd(LTR_nor_bySample)
LTR_can_sd <- sd(LTR_can_bySample)

LINE_normal_refseq <- chim_data[chim_data$repClass=="LINE" & chim_data$RefID!="." & chim_data$LIBRARY%%2==1,]
LINE_cancer_refseq <- chim_data[chim_data$repClass=="LINE" & chim_data$RefID!='.' & chim_data$LIBRARY%%2==0,]
SINE_normal_refseq <- chim_data[chim_data$repClass=="SINE" & chim_data$RefID!="." & chim_data$LIBRARY%%2==1,]
SINE_cancer_refseq <- chim_data[chim_data$repClass=="SINE" & chim_data$RefID!='.' & chim_data$LIBRARY%%2==0,]
LTR_normal_refseq <- chim_data[chim_data$repClass=="LTR" & chim_data$RefID!="." & chim_data$LIBRARY%%2==1,]
LTR_cancer_refseq <- chim_data[chim_data$repClass=="LTR" & chim_data$RefID!='.' & chim_data$LIBRARY%%2==0,]

LINE_nor_ref_bySample <- table(LINE_normal_refseq$LIBRARY) / LINE_ratio
LINE_can_ref_bySample <- table(LINE_cancer_refseq$LIBRARY) / LINE_ratio
SINE_nor_ref_bySample <- table(SINE_normal_refseq$LIBRARY) / SINE_ratio
SINE_can_ref_bySample <- table(SINE_cancer_refseq$LIBRARY) / SINE_ratio
LTR_nor_ref_bySample <- table(LTR_normal_refseq$LIBRARY) / LTR_ratio
LTR_can_ref_bySample <- table(LTR_cancer_refseq$LIBRARY) / LTR_ratio

LINE_nor_ref_mean <- mean(LINE_nor_ref_bySample)
LINE_can_ref_mean <- mean(LINE_can_ref_bySample)
SINE_nor_ref_mean <- mean(SINE_nor_ref_bySample)
SINE_can_ref_mean <- mean(SINE_can_ref_bySample)
LTR_nor_ref_mean <- mean(LTR_nor_ref_bySample)
LTR_can_ref_mean <- mean(LTR_can_ref_bySample)
LINE_nor_ref_sd <- sd(LINE_nor_ref_bySample)
LINE_can_ref_sd <- sd(LINE_can_ref_bySample)
SINE_nor_ref_sd <- sd(SINE_nor_ref_bySample)
SINE_can_ref_sd <- sd(SINE_can_ref_bySample)
LTR_nor_ref_sd <- sd(LTR_nor_ref_bySample)
LTR_can_ref_sd <- sd(LTR_can_ref_bySample)

LTR_normal_ref_sense <- chim_data[chim_data$repClass=="LTR" & chim_data$RefID!="." & chim_data$assXref=="s" & chim_data$LIBRARY%%2==1,]
LTR_cancer_ref_sense <- chim_data[chim_data$repClass=="LTR" & chim_data$RefID!='.' & chim_data$assXref=="s" & chim_data$LIBRARY%%2==0,]
LTR_normal_ref_anti <- chim_data[chim_data$repClass=="LTR" & chim_data$RefID!="." & chim_data$assXref=="as" & chim_data$LIBRARY%%2==1,]
LTR_cancer_ref_anti <- chim_data[chim_data$repClass=="LTR" & chim_data$RefID!='.' & chim_data$assXref=="as" & chim_data$LIBRARY%%2==0,]

LTR_nor_ref_sense_bySample <- table(LTR_normal_ref_sense$LIBRARY)
LTR_can_ref_sense_bySample <- table(LTR_cancer_ref_sense$LIBRARY)
LTR_nor_ref_anti_bySample <- table(LTR_normal_ref_anti$LIBRARY)
LTR_can_ref_anti_bySample <- table(LTR_cancer_ref_anti$LIBRARY)

LTR_nor_ref_sense_mean <- mean(LTR_nor_ref_sense_bySample)
LTR_can_ref_sense_mean <- mean(LTR_can_ref_sense_bySample)
LTR_nor_ref_anti_mean <- mean(LTR_nor_ref_anti_bySample)
LTR_can_ref_anti_mean <- mean(LTR_can_ref_anti_bySample)
LTR_nor_ref_sense_sd <- sd(LTR_nor_ref_sense_bySample)
LTR_can_ref_sense_sd <- sd(LTR_can_ref_sense_bySample)
LTR_nor_ref_anti_sd <- sd(LTR_nor_ref_anti_bySample)
LTR_can_ref_anti_sd <- sd(LTR_can_ref_anti_bySample)


#==============================================
# plotting
#==============================================


#error_bar.r
errorbars <- function (values, mid_points, errors) {
  if (length(mid_points) != length(errors)) {
    #    print "Number of mid points is different from number of errors!\n";
    return;
  }
  
  # calculate the width of error bars
  err_width <- (mid_points[2] - mid_points[1]) * 0.8 / 2;
  
  for (i in 1:length(mid_points)) {
    x_left <- mid_points[i] - err_width / 2.0;
    x_mid <- mid_points[i];
    x_right <- mid_points[i] + err_width / 2.0;
    y_top <- values[i] + errors[i];
    y_bottom <- values[i] - errors[i];
    segments(x_left, y_top, x_right, y_top);
    segments(x_left, y_bottom, x_right, y_bottom);
    segments(x_mid, y_top, x_mid, y_bottom);
  }
}

# plot for all chimeric transcripts
bar_heights <- matrix(c(LINE_nor_mean, LINE_can_mean, SINE_nor_mean, SINE_can_mean, LTR_nor_mean, LTR_can_mean), nrow=2)
colorset <- c("white", "darkblue");
par(mgp=c(4,1,0), mar=c(7,7,4,2)+0.1); # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colnames(bar_heights) <- c("LINE", "SINE", "LTR"); # attach bin names to data
rownames(bar_heights) <- c("Normal", "Cancer"); # attach group names to data
y_lim <- 200;
# making the plot
par(mfrow=c(1,2));
barplot(bar_heights, beside=T, col=colorset, axis.lty=1, ylim=c(0,y_lim), main="All Chimeric", xlab="Repeat Class", ylab="Normalized Chimeric Transcripts", cex.axis=0.9, cex.lab=1.2);
legend(2, y_lim, c("Normal", "Cancer"), colorset);
# plot error bars
mid_points <- c(1, 2, 4, 5, 7, 8) + 0.5;
SDs <- c(LINE_nor_sd, LINE_can_sd, SINE_nor_sd, SINE_can_sd, LTR_nor_sd, LTR_can_sd);
errorbars(bar_heights, mid_points, SDs);

# plot for chimeric transcripts associated with refSeq genes
bar_heights <- matrix(c(LINE_nor_ref_mean, LINE_can_ref_mean, SINE_nor_ref_mean, SINE_can_ref_mean, LTR_nor_ref_mean, LTR_can_ref_mean), nrow=2)
colnames(bar_heights) <- c("LINE", "SINE", "LTR"); # attach bin names to data
rownames(bar_heights) <- c("Normal", "Cancer"); # attach group names to data
y_lim <- 100;
barplot(bar_heights, beside=T, col=colorset, axis.lty=1, ylim=c(0,y_lim), main="RefSeq Only", xlab="Repeat Class", ylab="Normalized Chimeric Transcripts", cex.axis=0.9, cex.lab=1.2);
legend(2, y_lim, c("Normal", "Cancer"), colorset);
SDs <- c(LINE_nor_ref_sd, LINE_can_ref_sd, SINE_nor_ref_sd, SINE_can_ref_sd, LTR_nor_ref_sd, LTR_can_ref_sd);
errorbars(bar_heights, mid_points, SDs);

# plot for LTR chimeric transcripts by orientation
bar_heights <- matrix(c(LTR_nor_ref_anti_mean, LTR_nor_ref_sense_mean, LTR_can_ref_anti_mean, LTR_can_ref_sense_mean), nrow=2)
colnames(bar_heights) <- c("Normal", "Cancer"); # attach bin names to data
rownames(bar_heights) <- c("Antisense", "Sense"); # attach group names to data
y_lim <- 100;
par(mfrow=c(1,1));
colorset <- c("white", "green");
barplot(bar_heights, beside=T, col=colorset, axis.lty=1, ylim=c(0,y_lim), main="LTR Chimeric", xlab="Sample Type", ylab="Chimeric Transcripts", cex.axis=0.9, cex.lab=1.2);
legend(2, y_lim, c("Antisense", "Sense"), colorset);
SDs <- c(LTR_nor_ref_anti_sd, LTR_nor_ref_sense_sd, LTR_can_ref_anti_sd, LTR_can_ref_sense_sd);
errorbars(bar_heights, mid_points, SDs);

#================================================================================

infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/repFamilies_LINE.txt"
repFamilies_LINE <- read.table(infile, header=TRUE, sep="\t")
infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/repFamilies_SINE.txt"
repFamilies_SINE <- read.table(infile, header=TRUE, sep="\t")
infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/repFamilies_LTR.txt"
repFamilies_LTR <- read.table(infile, header=TRUE, sep="\t")

repFamilies <- c("L1", "L2", "Other")
repCount <- c(1001410, 474561, 68182+15605+9120+1113+532)
total_families <- length(repFamilies)
means <- numeric();
SDs <- numeric();
for (i in 1:total_families) {
  normal <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==1,]
  cancer <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==0,]
  nor_bySample <- table(normal$LIBRARY) / repCount[i]
  null_samples <- 63 - length(nor_bySample)
  nor_bySample <- c(nor_bySample, rep(0, null_samples))
  can_bySample <- table(cancer$LIBRARY) / repCount[i]
  null_samples <- 63 - length(can_bySample)
  can_bySample <- c(can_bySample, rep(0, null_samples))
  nor_mean <- mean(nor_bySample)
  can_mean <- mean(can_bySample)
  nor_sd <- sd(nor_bySample)
  can_sd <- sd(can_bySample)  
  means <- c(means, nor_mean, can_mean)
  SDs <- c(SDs, nor_sd, can_sd)
}
means[is.nan(means)]<-0
SDs[is.na(SDs)]<-0
bar_heights <- matrix(means, nrow=2)
colorset <- c("white", "darkblue");
par(mgp=c(4,1,0), mar=c(7,7,4,2)+0.1); # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colnames(bar_heights) <- repFamilies; # attach bin names to data
rownames(bar_heights) <- c("Normal", "Cancer"); # attach group names to data
# making the plot
par(mfrow=c(1,1));
y_lim <- max(means) * 1.5
barplot(bar_heights, beside=T, col=colorset, axis.lty=1, ylim=c(0,y_lim), main="LINE", xlab="Repeat Family", ylab="Frequency of Chimeric", cex.axis=0.9, cex.lab=1.2);
legend(2, y_lim, c("Normal", "Cancer"), colorset);
# plot error bars
mid_points <- c(1, 2, 4, 5, 7, 8) + 0.5;
errorbars(bar_heights, mid_points, SDs);

#====
# SINEs

repFamilies <- c("Alu", "MIR", "Other")
repCount <- c(1238897, 602609, 636+2550+5589+4295)
total_families <- length(repFamilies)
means <- numeric();
SDs <- numeric();
for (i in 1:total_families) {
  normal <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==1,]
  cancer <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==0,]
  nor_bySample <- table(normal$LIBRARY) / repCount[i]
  null_samples <- 63 - length(nor_bySample)
  nor_bySample <- c(nor_bySample, rep(0, null_samples))
  can_bySample <- table(cancer$LIBRARY) / repCount[i]
  null_samples <- 63 - length(can_bySample)
  can_bySample <- c(can_bySample, rep(0, null_samples))
  nor_mean <- mean(nor_bySample)
  can_mean <- mean(can_bySample)
  nor_sd <- sd(nor_bySample)
  can_sd <- sd(can_bySample)  
  means <- c(means, nor_mean, can_mean)
  SDs <- c(SDs, nor_sd, can_sd)
}
means[is.nan(means)]<-0
SDs[is.na(SDs)]<-0
bar_heights <- matrix(means, nrow=2)
colorset <- c("white", "darkblue");
par(mgp=c(4,1,0), mar=c(7,7,4,2)+0.1); # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colnames(bar_heights) <- repFamilies; # attach bin names to data
rownames(bar_heights) <- c("Normal", "Cancer"); # attach group names to data
# making the plot
par(mfrow=c(1,1));
y_lim <- max(means) * 1.5
barplot(bar_heights, beside=T, col=colorset, axis.lty=1, ylim=c(0,y_lim), main="SINE", xlab="Repeat Family", ylab="Frequency of Chimeric", cex.axis=0.9, cex.lab=1.2);
legend(2, y_lim, c("Normal", "Cancer"), colorset);
# plot error bars
mid_points <- c(1, 2, 4, 5, 7, 8) + 0.5;
errorbars(bar_heights, mid_points, SDs);


#====
# LTRs

repFamilies <- as.vector(repFamilies_LTR$repFamily)
repCount <- as.vector(repFamilies_LTR$count)
total_families <- length(repFamilies)
means <- numeric();
SDs <- numeric();
for (i in 1:total_families) {
  normal <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==1,]
  cancer <- chim_data[chim_data$repFamily==repFamilies[i] & chim_data$LIBRARY%%2==0,]
  nor_bySample <- table(normal$LIBRARY) / repCount[i]
  null_samples <- 63 - length(nor_bySample)
  nor_bySample <- c(nor_bySample, rep(0, null_samples))
  can_bySample <- table(cancer$LIBRARY) / repCount[i]
  null_samples <- 63 - length(can_bySample)
  can_bySample <- c(can_bySample, rep(0, null_samples))
  nor_mean <- mean(nor_bySample)
  can_mean <- mean(can_bySample)
  nor_sd <- sd(nor_bySample)
  can_sd <- sd(can_bySample)  
  means <- c(means, nor_mean, can_mean)
  SDs <- c(SDs, nor_sd, can_sd)
}
means[is.nan(means)]<-0
SDs[is.na(SDs)]<-0
bar_heights <- matrix(means, nrow=2)
colorset <- c("white", "darkblue");
par(mgp=c(4,1,0), mar=c(7,7,4,2)+0.1); # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colnames(bar_heights) <- repFamilies; # attach bin names to data
rownames(bar_heights) <- c("Normal", "Cancer"); # attach group names to data
# making the plot
par(mfrow=c(1,1));
y_lim <- max(means) * 1.5
barplot(bar_heights, beside=T, col=colorset, axis.lty=1, ylim=c(0,y_lim), main="LTR", xlab="Repeat Family", ylab="Frequency of Chimeric", cex.axis=0.9, cex.lab=1.2);
legend(2, y_lim, c("Normal", "Cancer"), colorset);
# plot error bars
mid_points <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17) + 0.5;
errorbars(bar_heights, mid_points, SDs);


