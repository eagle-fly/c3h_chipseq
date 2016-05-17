# check odds ratio of each mutation-chimera pair using Fisher's exact test
# Ying Zhang
# May 8, 2015

#threshold <- "relaxed"
threshold <- "stringent"

library(gplots)
setwd(paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, sep=""))

infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/mutations/muta_matrix.txt"
muta_matrix <- read.table(infile, header=TRUE, sep="\t")
infile <- "chim_matrix.txt"
chim_matrix <- read.table(infile, header=TRUE, sep="\t")

samples <- colnames(chim_matrix)
chimeras <- rownames(chim_matrix)
mutations <- rownames(muta_matrix)

infile <- "outliers.txt"
outlier_data <- read.table(infile, header=TRUE, sep="\t")
split <- unlist(strsplit(as.character(outlier_data$mc_pair), "-"))
muta <- split[seq(1, length(split) - 1, 2)]
chim <- split[seq(2, length(split), 2)]

outfile <- paste("../../results/", threshold, "/outlier_mutations.txt", sep="")
write.table(muta, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

# build presence matrix for outliers
presence <- numeric()
chim_unique <- unique(chim)
chim_samples <- character()
muta_samples <- character()

sideColorPalette <- rep(c("lightgray", "darkgray"), ceiling(length(chim_unique)/2))
rowSideColors <- character()
for (i in 1:length(chim_unique)) {
  ch <- chim_unique[i]
  # get the chimera index in chim_matrix
  ch_index <- which(chimeras==ch)
  presence <- rbind(presence, chim_matrix[ch_index,])

  # get the mutations paired with the current chimera
  mu_outliers <- muta[which(chim==ch)]
  mu_unique <- unique(mu_outliers)
  mu_index <- which(mutations %in% mu_unique)
  presence <- rbind(presence, muta_matrix[mu_index,])

  # output all mutations associated with the current chimera
  outfile <- paste("../../results/", threshold, "/outlier_muta_", ch, ".txt", sep="")
  write.table(mu_unique, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  
  # set the color of side bar
  rowSideColors <- c(rowSideColors, rep(sideColorPalette[i], length(mu_unique)+1))
}

# draw heatmap
x <- as.matrix(presence)
png(paste("../../results/", threshold, "/heatmap.png", sep=""), width=15, height=15, units="in", res=300)
heatmap.2(x, key=FALSE, margin=c(6,9),lmat=rbind(c(5,0,4), c(3,1,2)), lwid=c(1,1,10), lhei=c(1,10), 
          cexRow=1.2, cexCol=1.2, col=c("blue", "yellow"), RowSideColors=rowSideColors, trace="none", Rowv=NULL, Colv=NULL, dendrogram="none")
dev.off()

outfile <- paste("../../results/", threshold, "/outliers_presence.txt", sep="")
write.table(presence, outfile, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
