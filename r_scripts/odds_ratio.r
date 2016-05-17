# check odds ratio of each mutation-chimera pair using Fisher's exact test
# Ying Zhang
# May 8, 2015

threshold <- "relaxed"

infile <- "/Users/ying/Work/BCCRC/Projects/colon/data/mutations/muta_matrix.txt"
muta_matrix <- read.table(infile, header=TRUE, sep="\t")
infile <- paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, "/chim_matrix.txt", sep="")
chim_matrix <- read.table(infile, header=TRUE, sep="\t")

samples <- colnames(chim_matrix)
chimeras <- rownames(chim_matrix)
mutations <- rownames(muta_matrix)

odds_ratio <- matrix(nrow=length(mutations), ncol=length(chimeras))
p_values <- matrix(nrow=length(mutations), ncol=length(chimeras))
for (i in 1:length(mutations)) {
  for (j in 1:length(chimeras)) {
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
    odds_ratio[i,j] <- f$estimate
    p_values[i,j] <- f$p.value
  }
}

adjusted_p <- matrix(p.adjust(c(p_values), method="BH"), dim(p_values)[1], dim(p_values)[2])
