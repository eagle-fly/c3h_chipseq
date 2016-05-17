# check odds ratio of each mutation-chimera pair using Fisher's exact test
# Ying Zhang
# May 8, 2015

args<-commandArgs(TRUE)
setwd("/scratch/magerlab/yzhang/colon/hpc_relaxed")


infile <- "muta_matrix.txt"
muta_matrix <- read.table(infile, header=TRUE, sep="\t")
infile <- "chim_matrix.txt"
chim_matrix <- read.table(infile, header=TRUE, sep="\t")

samples <- colnames(chim_matrix)
chimeras <- rownames(chim_matrix)
mutations <- rownames(muta_matrix)

odds_ratio <- numeric()
p_values <- numeric()
j <- as.numeric(args[1])
for (i in 1:length(mutations)) {
muta <- muta_matrix[i,]
    chim <- chim_matrix[j,]
    subtraction <- muta - chim
    mu1_ch0 <- length(which(subtraction == 1))
    mu0_ch1 <- length(which(subtraction == -1))
    summation <- muta + chim
    mu1_ch1 <- length(which(summation == 2))
    mu0_ch0 <- length(which(summation == 0))
    mc_matrix <- matrix(c(mu1_ch1, mu0_ch1, mu1_ch0, mu0_ch0), 2, 2)
    f <- fisher.test(mc_matrix, alternative="greater")
    odds_ratio <- c(odds_ratio, as.numeric(f$estimate))
    p_values <- c(p_values, f$p.value)
}
result <- cbind(mutations, rep(chimeras[j], length(mutations)), odds_ratio, p_values)
colnames(result) <- c("mutations", "chimeras", "odds_ratio", "p_values")
outfile <- paste("pValues_chim", j, ".txt", sep="")
write.table(result, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

