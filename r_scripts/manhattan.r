threshold <- "stringent"
#threshold <- "relaxed"

setwd(paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, sep=""))
p_values <- numeric()
odds_ratio <- numeric()
chimeras <- character()
p_chim <- character()
mc_pair <- character()

infile <- "cancer_chim.txt"
chim <- read.table(infile, header=TRUE, sep="\t")
attach(chim)
chim_names <- paste(repClass, repFamily, repName, sep=":")
chim <- cbind(chim, chim_names)
detach(chim)
total_chim <- dim(chim)[1]

for (j in 1:total_chim) {
  infile <- paste("../../hpc_", threshold, "/pValues_chim", j, ".txt", sep="")
  data <- read.table(infile, header=TRUE, sep="\t")
  total_mutation <- dim(data)[1]
  repeatName <- as.character(chim[chim$repeatID == as.character(data$chimeras[1]), 9])
  chimeras <- c(chimeras, as.character(data$chimeras[1]))
  p_chim <- c(p_chim, repeatName)
  p_values <- c(p_values, data$p_values)
  odds_ratio <- c(odds_ratio, data$odds_ratio)
  mc_pair <- c(mc_pair, paste(as.character(data$mutations), as.character(data$chimeras), sep="-"))
}

#adjusted_p <- p.adjust(p_values, method="BH") # apply multiple test corrections

# prepare data for QQMAN
SNP <- mc_pair
repeatName <- p_chim
index_sorted <- order(repeatName)
index_sorted <- data.frame(index_sorted)
chim_index <- order(index_sorted$index_sorted)

CHR <- rep(chim_index, each=total_mutation)
BP <- rep(1, length(mc_pair))
P <- p_values
results <- data.frame(SNP, CHR, BP, P, odds_ratio)

# Manhattan plot
library(qqman)
png(paste("../../results/", threshold, "/manhattan.png", sep=""), width=10, height=6, units="in", res=300)
par(mar=c(7,7,5,2), mgp=c(5,1,0))
manhattan(results, xlab="Type of Chimera", suggestiveline= -log10(5e-02), genomewideline = -log10(1e-03), xaxt="n")
labels = p_chim[order(p_chim)]
axis(1, at=seq_along(labels), labels=FALSE)
text(x=seq_along(labels), y=par("usr")[3]-0.2, srt=45, adj = 1, labels = labels, cex=0.5, xpd = TRUE)
labels = chimeras[order(p_chim)]
axis(3, at=seq_along(labels), labels=FALSE)
text(x=seq_along(labels), y=par("usr")[4]+0.2, srt=-45, adj = 1, labels = labels, cex=0.5, xpd = TRUE)
dev.off()

outliers <- results[results$P <= 0.001, c(1,4,5)]
colnames(outliers) <- c("mc_pair", "p_value", "odds_ratio")
split <- unlist(strsplit(as.character(outliers$mc_pair), "-"))
mutation <- split[seq(1, length(split) - 1, 2)]
repeatID <- split[seq(2, length(split), 2)]
outliers <- cbind(outliers, mutation, repeatID)
outliers <- merge(outliers, chim)
outliers <- outliers[order(outliers$p_value),]
outfile <- paste("../../results/", threshold, "/outliers.txt", sep="")
write.table(outliers, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
