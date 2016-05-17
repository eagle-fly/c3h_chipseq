# Compare the library size (mapped reads) between cancer and normal samples
# Ying Zhang
# June 5, 2015

infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/lib_sizes.txt"
lib_sizes <- read.table(infile, header=TRUE, sep="\t")

normal <- lib_sizes[lib_sizes$Sample%%2==1, 2]
cancer <- lib_sizes[lib_sizes$Sample%%2==0, 2]
sample_list <- lib_sizes$Sample

# boxplot of overall sample size difference
y_lim <- max(c(normal, cancer))
png(paste("../../results/libSize_boxplot.png", sep=""), width=6, height=6, units="in", res=300)
par(mgp=c(4,1,0), mar=c(6,6,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("orange", "blue")
boxplot(list(normal, cancer), range=1, col=colorset,
        ylab="Mapped reads in sample", 
        cex.axis=0.9, cex.lab=1.2, names=c("Normal", "Cancer"))
dev.off()


# barplot of library size between normal and cancer

png(paste("../../results/libSize_barplot.png", sep=""), width=16, height=8, units="in", res=300)
par(mgp=c(3,1,0), mar=c(5,5,4,2)+0.1) # mgp sets the distance of axis labels to the axis; mar sets the size of figure margin
colorset <- c("darkblue", "orange")
barplot(rbind(cancer, normal), beside=TRUE, col=colorset, xlab="Samples", ylab="Sample Size")
axis(1, tick=FALSE, labels=FALSE)
text(x=sort(c(seq(1,187,3), seq(2,188,3)))+0.3, y=par("usr")[3]-0.2, srt=90, adj = 0.8, labels = sample_list, cex=0.5, xpd = TRUE)
legend(1, y_lim, c("Cancer", "Normal"), colorset)
dev.off()

# load chimera data

#threshold <- "stringent"
threshold <- "relaxed"
te <- "SINE"

setwd(paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, sep=""))
infile <- "sorted_chim/chim_all.txt"
chim_data <- read.table(infile, header=TRUE, sep="\t")
load("/Users/ying/Work/BCCRC/Projects/colon/artem_scripts/RM/RepeatMasker.Rdata")


TE_total <- 1407573768
LINE_total <- claStats[3,3]
SINE_total <- claStats[15,3]
LTR_total <- claStats[6,3]

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

LINE_chimeras <- as.numeric(LINE_nor_bySample)
SINE_chimeras <- as.numeric(SINE_nor_bySample)
LTR_chimeras <- as.numeric(LTR_nor_bySample)

if (te == "LINE") {
  chimeras <- LINE_chimeras
} else if (te == "SINE") {
  chimeras <- SINE_chimeras
} else if (te == "LTR") {
  chimeras <- LTR_chimeras
}
r <- cor(normal, chimeras)
png(paste("../../results/", threshold, "/libSize_chim_", te, ".png", sep=""), width=6, height=6, units="in", res=300)
plot(normal, chimeras, xlab="Mapped Reads in Sample", ylab= paste(te, " Chimeras in Sample", sep=""),
     main=paste("Normal (r = ", format(r, digits=2), ")", sep=""))
dev.off()
