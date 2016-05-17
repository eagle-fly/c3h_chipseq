# Build the presence-absence matrices for Mutations vs. Samples and Chimera vs. Samples
# Ying Zhang
# April 27, 2015

threshold <- "relaxed"

setwd(paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, sep=""))
infile <- "cancer_chimeras.txt"
chim_data <- read.table(infile, header=TRUE, sep="\t")
infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/mutations/mutations.txt"
muta_data <- read.table(infile, header=TRUE, sep="\t")
infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/list.txt"
samples <- read.table(infile, header=FALSE, sep="\t")
samples <- samples[,1]
samples <- samples[samples %% 2 == 0] # take only cancer samples

mutations <- as.vector(unique(muta_data$Gene))
chimeras <- as.vector(unique(chim_data$repeatID))
outfile <- "unique_mutations.txt"
write.table(mutations, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
outfile <- "unique_chimeras.txt"
write.table(chimeras, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

# build the mutation presence matrix
muta_matrix <- matrix(rep(0, length(mutations) * length(samples)), nrow=length(mutations),ncol=length(samples))
colnames(muta_matrix) <- samples
rownames(muta_matrix) <- mutations
for (j in 1:length(samples)) {
  genes <- unique(as.vector(muta_data[muta_data$Sample == samples[j], 2]))
  gene_index <- which(mutations %in% genes)
  muta_matrix[gene_index, j] <- 1
}

# build the chimera presence matrix
chim_matrix <- matrix(rep(0, length(chimeras) * length(samples)), nrow=length(chimeras),ncol=length(samples))
colnames(chim_matrix) <- samples
rownames(chim_matrix) <- chimeras
for (j in 1:length(samples)) {
  TEs <- unique(as.vector(chim_data[chim_data$Sample == samples[j], 2]))
  TE_index <- which(chimeras %in% TEs)
  chim_matrix[TE_index, j] <- 1
}

outfile <- "muta_matrix.txt"
write.table(muta_matrix, outfile, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
outfile <- "chim_matrix.txt"
write.table(chim_matrix, outfile, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)


