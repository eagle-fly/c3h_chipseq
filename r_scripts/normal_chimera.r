# check if a given chimeric transcript is cancer-specific based on Fisher's exact test
# Ying Zhang
# May 6, 2015


threshold <- "stringent"

setwd(paste("/Users/ying/Work/BCCRC/Projects/colon/data/", threshold, sep=""))
infile <- "chimeras.txt"
chim_data <- read.table(infile, header=TRUE, sep="\t")

chim_TEs <- unique(as.vector(chim_data$repeatID))
odds_ratio <- numeric()
p_values <- numeric()
for (te in chim_TEs) {
  samples <- unique(as.vector(chim_data[chim_data$repeatID == te, 1]))
  cancer <- length(which(samples %% 2 == 0))
  normal <- length(which(samples %% 2 == 1))
  cancer_null <- 63 - cancer
  normal_null <- 63 - normal
  c <- matrix(c(normal, cancer, normal_null, cancer_null), 2, 2)
  f <- fisher.test(c, alternative="greater")
  odds_ratio <- c(odds_ratio, f$estimate)
  p_values <- c(p_values, f$p.value)
}
adjusted_p <- p.adjust(p_values, method="BH")
results <- data.frame(chim_TEs, odds_ratio, p_values, adjusted_p)
outfile <- "normal_chimeras_stat.txt"
write.table(results, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

normal_chimTEs <- as.vector(results[results$adjusted_p <= 0.05, 1])
normal_chim <- chim_data[chim_data$repeatID %in% normal_chimTEs,]
outfile <- "normal_chimeras.txt"
write.table(normal_chim, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

# output refIDs
refIDs <- as.character(unique(normal_chim$refID))
refIDs <- refIDs[refIDs != "."]
outfile <- "normal_chimeras_refID.txt"
write.table(refIDs, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

